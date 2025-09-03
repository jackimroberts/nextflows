#!/usr/bin/env nextflow

// Include subworkflows and shared help
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { getSharedHelp } from '../modules/shared_help'
include { handleWorkflowCompletion } from '../shared_bin/workflow_completion_handler.nf'

params.miniaturize = false
params.fastq_source = true
params.sample_table = "*.txt"

/*
 * Things I'd like to add:
 ** qc throughout
 ** create containerized version
 ** deploy on cloud
*/

params.help = false
if (params.help) {
        println """
        ATAC-seq nextflow pipeline
${getSharedHelp()}
        """
        exit 0
}

workflow {

	// FASTQ preprocessing (sample table creation, download, decompress, miniaturize)
	FASTQ_PREPROCESSING(
		Channel.fromPath(params.sample_table),
		params.fastq_source, 
		params.miniaturize
	)
	
	sample_sheet = FASTQ_PREPROCESSING.out.sample_sheet
	fastq_list = FASTQ_PREPROCESSING.out.fastq_files
		
	// pair up fastq files by sample id and run, keeping runs separate
 	sample_sheet
     		| combine(fastq_list)
      		| filter { sample_id, sample_name, condition, fastq_file ->
          		fastq_file.name.startsWith(sample_id + "_")
      		}
      		| map { sample_id, sample_name, condition, fastq_file ->
          		def run_part = fastq_file.name.replaceAll("^${sample_id}_", "").replaceAll("[-_]?[LR][12].*", "")
          		def meta = [id: sample_id, name: sample_name, condition: condition, run: run_part]
          		[[meta.id, meta.run], meta, fastq_file]
      		}
      		| groupTuple()  // Group by the [id, run] tuple
      		| map { id_run_tuple, meta_list, fastq_files ->
          		[meta_list[0], fastq_files]
      		}
	// Trim and align individual sequencing runs
		| adapter_trim 
		| bowtie_align
	// Combine samples by sample id
		| map { meta, bam_file ->
			[meta.id, meta, bam_file]
		}
		| groupTuple
	// Merge bam files, filter, and calculate depth
		| map { sample_id, meta_list, bam_files ->
			[meta_list[0], bam_files]
	 	}
		| merge_run_bams
		| filter_bams
		| set { filtered_bams}

	// Call peaks from each sample
	filtered_bams
		| call_peaks
		| set { peak_outputs }
	// Find consensus peaks
	peak_outputs.peak_files
		| collect
		| merge_peaks
	// Obtain counts from consensus peaks
		| combine(peak_outputs.coverage_files)
		| count_under_peaks
	// DE analysis
	// group by condition

	// Calculate read depth for each bam
	filtered_bams
		| measure_depth
		| map { meta, bam_file, depth ->
          		def depth_value = depth.trim() as Integer
   
          		[meta + [depth: depth_value], bam_file]
		}
		| set { filtered_bams_with_depth }
    	// Find minimum depth across all samples for scaling
  	filtered_bams_with_depth
      		| map { meta, bam_file -> meta.depth }
      		| min()
  	// Combine for bigwig creation
      		| combine(filtered_bams_with_depth)
      		| bam_to_bigwig
	
	// run Tim's pipeline
	filtered_bams_with_depth
		| map { meta, bam_file ->
			[meta.condition, meta, bam_file]
		}
		| groupTuple
		| map { condition, meta_list, bam_files ->
			[meta_list[0], bam_files]
		}
		| create_multimacs_run
		| collect
		| run_multimacs
}

workflow.onComplete {
    handleWorkflowCompletion()
}


/*
 * adapter trimming
 */
process adapter_trim {
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path('*fq')
	script:
		"""
          	#!/bin/bash
          	# TruSeq adapters for Active Motif ATAC-seq
		adaptf=CTGTCTCTTATACACATCT
		adaptr=CTGTCTCTTATACACATCT

		# load cutadapt module
		module load cutadapt

		#cutadapt
		echo "=== adapter trimming"
		which cutadapt
		cutadapt -O 1 --nextseq-trim=20 -m 1 -a \$adaptf -A \$adaptr \\
			-o ${meta.id}_${meta.run}.1.fq -p ${meta.id}_${meta.run}.2.fq \\
			$fastqs

		"""
}

/*
 * alignment with bowtie2
 */
process bowtie_align {
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path('*raw.bam')
	script:
		def r1 = fastqs.find { it.name.contains('.1.fq') }
		def r2 = fastqs.find { it.name.contains('.2.fq') }
		"""
		#!/bin/bash

		# load required modules
		module load bowtie2
		module load samtools

		echo "=== align paired-end reads with bowtie2"
		which bowtie2
		which samtools

		bowtie2 --local --very-sensitive --no-mixed --no-discordant -p \$(nproc) \\
			-x ${params.genome_index} -1 ${r1} -2 ${r2} | \\
			samtools view -bS -o ${meta.id}_${meta.run}.raw.bam
		"""
}

/*
 * Combine BAM files from the same sample
 */
process merge_run_bams {
    input:
        tuple val(meta), path(bam_files)
    output:
        tuple val(meta), path("${meta.id}.merged.bam")
    script:
        """
        module load samtools
        
        if [ \$(echo $bam_files | wc -w) -gt 1 ]; then
            samtools merge ${meta.id}.merged.bam $bam_files
        else
            # Use readlink to resolve any symbolic links to actual file
            actual_file=\$(readlink -f $bam_files)
            ln -s "\$actual_file" ${meta.id}.merged.bam
        fi
        """
}

/*
 * filter bams for duplicates, unaligned, mitochondria, blacklist
 */
process filter_bams {
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path('*filtered.bam')
	script:
		"""
		sh ${projectDir}/bin/filter_bams.sh $bam ${meta.id} ${params.genome_blacklist}
		"""
}
/*
 * Call peaks from fragment ends
 * This strategy is specific to ATAC-seq data
 * Assume the fragment's 5' end is the open chromatin as it's the cut site
 * Take the 5' fragment ends as transcript locations with shift 100 ext 200
 */
process call_peaks {
	publishDir 'peaks', mode: 'copy'
	input:
		tuple val(meta), path(bam)
	output:
		path("*narrowPeak"), emit: peak_files
		tuple val(meta), path(bam), path("${meta.id}.bed"), emit: coverage_files 
	script:
		"""
		#!/bin/bash
		
		# load required modules
		module load bedtools
		module load macs
		bedtools --version
		macs2 --version
		
		# Convert to bed file with each read as an individual fragment
		echo "=== Convert bam to bed"
		bedtools bamtobed -i $bam > "${meta.name}.bed"
		
		# Call peaks on 5' end
		echo "=== Call peaks with macs2"
		macs2 callpeak -t "${meta.name}.bed" -f BED -n ${meta.name} \\
    			--keep-dup all --qvalue 0.01 \\
    			--shift -100 --extsize 200 \\
    			--min-length 100 --gsize 2500000000 2> "${meta.name}.macs.log"
		"""
}

/*
 * Merge peaks from multiple samples into consensus peak set
 * Strategy: Keep peaks found in at least 2 samples, merge nearby peaks
 */
process merge_peaks {
	publishDir 'peaks', mode: 'copy'
	input:
		path(peak_files)
	output:
		path "merged.bed"
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load bedtools
		bedtools --version

		echo "=== Creating consensus peak set from \$(echo $peak_files | wc -w) samples"
		
		# Step 1: Combine all peaks and sort by genomic coordinates
		cat $peak_files | sort -k1,1 -k2,2n > all_peaks.bed
		
		# Step 2: Find overlapping peaks across samples
		# -a: query intervals (all peaks)
		# -b: database intervals (all peak files) 
		# -c: count overlaps with full requirement (-f 1)
		# -f 1: require 100% overlap of query interval
		bedtools intersect -a all_peaks.bed -b $peak_files -c -f 1 > peaks_with_counts.bed
		
		# Step 3: Keep only peaks found in at least 2 samples
		# Column 11 contains the overlap count from bedtools intersect -c
		awk '\$11 >= 2' peaks_with_counts.bed | sort -k1,1 -k2,2n > filtered_peaks.bed
		
		# Step 4: Merge peaks within 100bp of each other
		bedtools merge -i filtered_peaks.bed -d 100 > merged.bed
		
		echo "Final consensus peak set: \$(wc -l < merged.bed) regions"
		
		# Clean up intermediate files
		rm all_peaks.bed peaks_with_counts.bed filtered_peaks.bed
		"""
}

/*
 * Count reads overlapping consensus peaks for each sample
 * Uses both bedtools and featureCounts approaches for comparison
 */
process count_under_peaks {
	publishDir 'peak_counts', mode: 'copy'
	input:
		tuple path(merged_peaks), val(meta), path(bam), path(bed)
	output:
		path "${meta.id}_counts_*.txt"
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load bedtools
		module load subread  # for featureCounts
		bedtools --version
		featureCounts -v

		echo "=== Counting reads under consensus peaks for sample ${meta.id}"
		
		# Method 1: Using bedtools intersect
		# Sort the sample bed file by genomic coordinates
		sort -k1,1 -k2,2n $bed > ${meta.id}_sorted.bed
		
		# Count overlaps between consensus peaks and sample fragments
		# -a: consensus peaks (query)
		# -b: sample fragments (database) 
		# -wa: write original entry from -a file
		# -c: count overlaps instead of writing them all
		bedtools intersect -a $merged_peaks -b ${meta.id}_sorted.bed -wa -c > ${meta.id}_counts_bedtools.txt
		
		# Method 2: Using featureCounts (more standard for RNA-seq/ChIP-seq)
		# Convert consensus peaks to SAF format for featureCounts
		# SAF format: GeneID, Chr, Start, End, Strand
		awk 'BEGIN{OFS="\\t"} {print \$1"_"\$2"_"\$3, \$1, \$2, \$3, "."}' $merged_peaks > consensus_peaks.saf
		
		# Count reads overlapping peaks using featureCounts
		# -a: annotation file (SAF format)
		# -o: output file
		# -p: count fragments (paired-end mode)
		featureCounts -p -a consensus_peaks.saf -o ${meta.id}_counts_featureCounts.txt $bam
		
		# Clean up intermediate files
		rm ${meta.id}_sorted.bed consensus_peaks.saf
		"""
}

/*
 * Calculate bam sizes for scaling
 * append sequencing depth to output
 */
process measure_depth {
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path(bam), stdout
	script:
		"""
		module load samtools
		
		## Calculate size for eventual scaling
		echo \$(samtools view -f 64 -c $bam) | tr -d '\n'

		"""
}

/*
 * convert bam files to bigwig
 */
process bam_to_bigwig {
	publishDir 'bigwigs', mode: 'copy'
	input:
		tuple val(min_depth), val(meta), path(bam)
	output:
		path "*bw"
	script:
		"""
		#!/bin/bash
		
		# load required modules - load deeptools first to avoid Python conflicts
		module purge
		module load deeptools
		module load samtools

		# create scale based on depth and minimum depth of all samples
		scale=\$(echo "scale=5; $min_depth/${meta.depth}" | bc)

		# index bam and generate bigwig
		echo "=== generate bigwig from bam, scaled to smallest bam in group"
		which samtools
		which deeptools

		samtools index $bam
		bamCoverage -b $bam --scaleFactor \$scale -o ${meta.name}.bw
		"""
}

/*
 * run filtered bam files through Tim's pipeline
 */
process create_multimacs_run {
	maxForks 1 // Each condition is added sequentially
	input: 
		tuple val(meta), val(bam_files)
	output: 
		stdout emit: "done"
	script:
		"""
		# Move run template to launchDir
		if ! test -f ${launchDir}/output/multimacs_command.sh; then
			mkdir -p ${launchDir}/output
  			cp ${projectDir}/bin/multimacs_run_template.txt ${launchDir}/output/multimacs_command.sh
		fi


		# Create text with the format:
			# --chip /file1,/file2 --name condition1
		clean_bams=\$(echo "${bam_files}" | tr -d '[]' | tr ', ' ',')
		chip_details="multirep_macs2_pipeline.pl --chip \$clean_bams --name ${meta.condition}"
		
		# Place run-specific details in the multimacs script
		sed -i "s~multirep_macs2_pipeline.pl~\$chip_details~" ${launchDir}/output/multimacs_command.sh
		echo "done"
		"""
}
process run_multimacs {
	//publishDir 'multimacs_pipeline', mode: 'move'
	input:
		val _
	output:
		path "*"
	script:
		"""
		${launchDir}/output/multimacs_command.sh
		"""
}



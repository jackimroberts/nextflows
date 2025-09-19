#!/usr/bin/env nextflow

params.miniaturize = false
params.fastq_source = true
params.sample_table = "*.txt"
params.outputDir = "output"

//ATAC-specific. Does not remove unmated. Removes, duplicate, unaligned, and secondary.
params.filter = [flag: 3332]

// Include subworkflows and shared help
include { getSharedHelp } from '../modules/shared_help'
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { MERGE_SEQUENCING_RUNS } from '../subworkflows/merge_sequencing_runs.nf'
include { CREATE_SCALED_BIGWIGS } from '../subworkflows/create_scaled_bigwigs.nf'
include { adapter_trim; filter_bams; collect_qc; multiqc } from '../modules/shared_processes.nf'
include { WorkflowCompletion } from '../subworkflows/workflow_complete.nf'

params.help = false
if (params.help) {
        println """
        ATAC-seq nextflow pipeline
${getSharedHelp()}

        ATAC-seq specific parameters:
        --macs2.keepDup, .qvalue, .shift, .extsize, .minLength, .gsize, 
        .args (additional command line arguments)
                MACS2 peak calling parameters (see macs2 callpeak --help)
                
        --adapters.forward, .reverse, .overlap, .nextseqTrim, .minLength, .args
                Adapter trimming parameters (see cutadapt --help)

        --filter.flag integer
		--filter.flag 3332 or --filter ["flag": 3332]
                SAM filtering flags (default: 3332 for ATAC-seq)
        """
        exit 0
}

workflow {

	// FASTQ preprocessing (sample table creation, download, decompress, miniaturize, fastq pairing/grouping)
	FASTQ_PREPROCESSING(
		Channel.fromPath(params.sample_table),
		params.fastq_source, 
		params.miniaturize
	)
	
	grouped_fastqs = FASTQ_PREPROCESSING.out.grouped_fastqs
	// Trim and align individual sequencing runs
		| adapter_trim 
		| bowtie_align
	
	// Merge samples across sequencing runs
	MERGE_SEQUENCING_RUNS(bowtie_align.out)
	
	// Filter merged BAM files
	merged_bams = MERGE_SEQUENCING_RUNS.out.merged_bams
		| filter_bams
		| map { meta, bam, bai -> [meta, bam] }
		| set { filtered_bams }

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
	// Collect QC files and generate MultiQC report
		| collect
		| collect_qc
		| multiqc
	// DE analysis
	// group by condition

	// Create scaled bigwig files for visualization
	CREATE_SCALED_BIGWIGS(filtered_bams)
	
	// run Tim's pipeline
	filtered_bams
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
	// Runs on success, cancel or fail
	WorkflowCompletion()
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

		echo "====== PROCESS_SUMMARY"
		echo "====== BOWTIE_ALIGN ======"
		echo "Strategy: Align paired-end reads with bowtie2"
		echo "bowtie2 version \$(bowtie2 --version | head -1 | awk '{print \$NF}')"
		echo "\$(samtools --version | head -1)"
		echo "Parameters: --local --very-sensitive --no-mixed --no-discordant"
		echo "Reference: ${params.genome_index}"
		echo "====== BOWTIE_ALIGN ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== aligning $r1/.2.fq" 
		bowtie2 --local --very-sensitive --no-mixed --no-discordant -p ${task.cpus} \\
			-x ${params.genome_index} -1 ${r1} -2 ${r2} \\
			2> ${meta.id}${meta.run}_bowtie.log | \\
			samtools view -bS -o ${meta.id}${meta.run}.raw.bam

		# Extract key stats from bowtie stderr (last 6 lines)
		echo "=== Alignment Results"
		tail -6 ${meta.id}${meta.run}_bowtie.log
		
		echo "=== Output Files"
		bam_count=\$(samtools view -c ${meta.id}${meta.run}.raw.bam)
		echo "${meta.id}${meta.run}.raw.bam: \$bam_count reads aligned"
		"""
}

/*
 * Call peaks from fragment ends
 * This strategy is specific to ATAC-seq data
 * Assume the fragment's 5' end is the open chromatin as it's the cut site
 * Take the 5' fragment ends as transcript locations with shift 100 ext 200
 */
process call_peaks {
	publishDir "${params.outputDir}/peaks", mode: 'copy', pattern: '*_peaks.narrowPeak'
	params.macs2 = [
		keepDup: 'all',              // keep duplicate reads
		qvalue: 0.01,                // q-value threshold
		shift: -100,                 // shift reads (ATAC-seq specific)
		extsize: 200,                // extend reads to this length
		minLength: 100,              // minimum peak length
		gsize: 2500000000,           // effective genome size
		args: ''                     // custom MACS2 arguments
	]
	input:
		tuple val(meta), path(bam)
	output:
		path("${meta.name}_peaks.narrowPeak"), emit: peak_files
		tuple val(meta), path(bam), path("${meta.name}.bed"), emit: coverage_files 
	script:
		"""
		#!/bin/bash
		
		# load required modules
		module load bedtools
		module load macs

		echo "====== PROCESS_SUMMARY"
		echo "====== CALL_PEAKS ======"
		echo "Strategy: Call peaks from fragment 5' ends"
		echo "Convert bam files to bed using bedtools \$(bedtools --version)"
		echo "Call peaks with macs2 \$(macs2 --version)"
		echo "Parameters:"
		echo "  keepDup: ${params.macs2.keepDup}"
		echo "  qvalue: ${params.macs2.qvalue}"
		echo "  shift: ${params.macs2.shift} (ATAC-seq specific)"
		echo "  extsize: ${params.macs2.extsize}"
		echo "  minLength: ${params.macs2.minLength}"
		echo "  gsize: ${params.macs2.gsize}"
		if [ -n "${params.macs2.args}" ]; then
			echo "  Custom args: ${params.macs2.args}"
		fi
		echo "====== CALL_PEAKS ======"
		echo "====== PROCESS_SUMMARY"
		
		# Convert to bed file with each read as an individual fragment
		echo "=== Convert $bam to ${meta.name}.bed"
		bedtools bamtobed -i $bam > "${meta.name}.bed"
		
		# Confirm bed file creation
		bed_lines=\$(wc -l < "${meta.name}.bed")
		bed_size=\$(du -h "${meta.name}.bed" | cut -f1)
		echo "Created ${meta.name}.bed: \$bed_lines fragments (\$bed_size)"
		
		# Call peaks on 5' end
		echo "=== Call peaks from ${meta.name}.bed"
		macs2 callpeak \\
			-t "${meta.name}.bed" \\
			-f BED \\
			-n ${meta.name} \\
			--keep-dup ${params.macs2.keepDup} \\
			--qvalue ${params.macs2.qvalue} \\
			--shift ${params.macs2.shift} \\
			--extsize ${params.macs2.extsize} \\
			--min-length ${params.macs2.minLength} \\
			--gsize ${params.macs2.gsize} \\
			${params.macs2.args} \\
			2> "${meta.name}.macs.log"
		
		# Extract key stats from MACS2 log
		grep "total tags in treatment" "${meta.name}.macs.log"
		tail -1 "${meta.name}.macs.log"
		
		# Ensure peak file exists (create empty if MACS2 found no peaks)
		if [ ! -f "${meta.name}_peaks.narrowPeak" ]; then
			touch "${meta.name}_peaks.narrowPeak"
		fi
		
		# Confirm peak file creation
		peak_count=\$(wc -l < "${meta.name}_peaks.narrowPeak")
		echo "Called \$peak_count peaks"
		"""
}

/*
 * Merge peaks from multiple samples into consensus peak set
 * Strategy: Keep peaks found in at least 2 samples, merge nearby peaks
 */
process merge_peaks {
	publishDir "${params.outputDir}/peaks", mode: 'copy'
	input:
		path(peak_files)
	output:
		path "merged.bed"
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load bedtools
		
		echo "====== PROCESS_SUMMARY"
		echo "====== MERGE_PEAKS ======"
		echo "Strategy: Create consensus peaks from multiple samples"
		echo "bedtools \$(bedtools --version)"
		echo "Keep peaks found in >=2 samples, merge nearby peaks (100bp)"
		echo "====== MERGE_PEAKS ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== Creating consensus peak set from \$(echo $peak_files | wc -w) samples"
		echo "${peak_files}"

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
		if [ -s peaks_with_counts.bed ]; then
			awk '\$11 >= 2' peaks_with_counts.bed | sort -k1,1 -k2,2n > filtered_peaks.bed
		else
			touch filtered_peaks.bed
		fi
		
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
	publishDir "${params.outputDir}/counts", mode: 'copy'
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
		
		echo "====== PROCESS_SUMMARY"
		echo "====== COUNT_UNDER_PEAKS ======"
		echo "Strategy: Count reads overlapping consensus peaks"
		echo "bedtools \$(bedtools --version)"
		echo "subread (featureCounts) \$(subread --version)"
		echo "Uses both bedtools intersect and featureCounts methods"
		echo "====== COUNT_UNDER_PEAKS ======"
		echo "====== PROCESS_SUMMARY"

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
		
		echo "=== Count Results ==="
		# Total counts from bedtools method (sum column 4)
		bedtools_total=\$(awk '{sum+=\$4} END {print sum}' ${meta.id}_counts_bedtools.txt)
		echo "bedtools intersect: \$bedtools_total total counts"
		
		# Total counts from featureCounts method (sum column 7, skip header lines)
		featureCounts_total=\$(tail -n +3 ${meta.id}_counts_featureCounts.txt | awk '{sum+=\$7} END {print sum}')
		echo "featureCounts: \$featureCounts_total total counts"
		
		# Clean up intermediate files
		rm ${meta.id}_sorted.bed consensus_peaks.saf
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
		echo "====== PROCESS_SUMMARY"
		echo "====== CREATE_MULTIMACS_RUN ======"
		echo "Strategy: Prepare multimacs command for peak calling"
		echo "Introduce --chip [bams] --name [condition] to template"
		echo "====== CREATE_MULTIMACS_RUN ======"
		echo "====== PROCESS_SUMMARY"
		
		# Move run template to launchDir
		if ! test -f ${launchDir}/output/multimacs_command.sh; then
			mkdir -p ${launchDir}/output
  			cp ${projectDir}/bin/multimacs_run_template.txt ${launchDir}/output/multimacs_command.sh
		fi


		# Create text with the format:
			# --chip /file1,/file2 --name condition1
		clean_bams=\$(echo "${bam_files}" | tr -d '[]' | sed 's/, /,/g')
		chip_details="multirep_macs2_pipeline.pl --chip \$clean_bams --name ${meta.condition}"
		
		# Place run-specific details in the multimacs script
		sed -i "s~multirep_macs2_pipeline.pl~\$chip_details~" ${launchDir}/output/multimacs_command.sh
		echo "\$chip_details"
		"""
}
process run_multimacs {
	//publishDir 'multimacs_pipeline', mode: 'move'
	errorStrategy 'ignore'
	input:
		val _
	output:
		path "*", optional: true
	script:
		"""
		module use /uufs/chpc.utah.edu/common/home/hcibcore/Modules/modulefiles

		module load multirepchipseq
		module load parallel
		module load bedtools

		echo "====== PROCESS_SUMMARY"
		echo "====== RUN_MULTIMACS ======"
		echo "Strategy: Execute multimacs pipeline for peak analysis"
		echo "multirepchipseq \$(multirepchipseq --version)"
		echo "parallel \$(parallel --version)"
		echo "bedtools \$(bedtools --version)"
		echo "run parameters are in launchDir/output"
		echo "====== RUN_MULTIMACS ======"
		echo "====== PROCESS_SUMMARY"

		${launchDir}/output/multimacs_command.sh ${task.cpus}
		"""
}



#!/usr/bin/env nextflow

// Include subworkflows and shared help
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { WORKFLOW_COMPLETION } from '../subworkflows/workflow_completion.nf'
include { getSharedHelp } from '../modules/shared_help'

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
		| map { meta, bam_file, depth ->
          		def depth_value = depth.trim().split('\n').last()
          		def more_meta = meta + [depth: depth_value]
          		[more_meta, bam_file]
		}
		| set { filtered_bams_with_depth }

  	// Find minimum depth across all samples for scaling
	// Create scaled bigwigs
  	filtered_bams_with_depth
      		| map { meta, bam_file -> meta.depth as Integer }
      		| min()
      		| combine(filtered_bams_with_depth) // [min_depth, meta, bam_file]
      		| bam_to_bigwig
		
	/*
	// callpeaks, merge
	
	// counts under peaks
	// DE analysis
	// group by condition
	bam_sample_sheet|
	| map { row -> [row.cond,row.file]}
	| group_tuple 
	| run_multimacs
	*/

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
	
	// Trigger completion after all processes finish
	run_multimacs.out
		| mix(bam_to_bigwig.out)
		| WORKFLOW_COMPLETION
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
 * also append sequencing depth to output
 */
process filter_bams {
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path('*filtered.bam'), stdout
	script:
		"""
		sh ${projectDir}/bin/filter_bams.sh $bam ${meta.id} ${params.genome_blacklist}
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



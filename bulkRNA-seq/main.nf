#!/usr/bin/env nextflow

// Include subworkflows and shared help
include { getSharedHelp } from '../modules/shared_help'
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { CREATE_SCALED_BIGWIGS } from '../subworkflows/create_scaled_bigwigs.nf'
include { MERGE_SEQUENCING_RUNS } from '../subworkflows/merge_sequencing_runs.nf'
include { adapter_trim; filter_bams } from '../modules/shared_processes.nf'
include { WorkflowCompletion } from '../subworkflows/workflow_complete.nf'

params.miniaturize = false
params.fastq_source = true
params.sample_table = "*.txt"

params.help = false
if (params.help) {
        println """
        ATAC-seq nextflow pipeline
${getSharedHelp()}
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
		| fastqc1
		| extract_umi
		| fastqc2
		| adapter_trim 
		| star_align
	
	// Merge samples across sequencing runs
	MERGE_SEQUENCING_RUNS(star_align.out)
	
	// Filter merged BAM files
	merged_bams = MERGE_SEQUENCING_RUNS.out.merged_bams
		| filter_bams
		| dedup_umi
		| dedup_qc
		| set { filtered_bams }

	// Count features in each sample
	filtered_bams
		| count_features
		| collect
		| run_DESeq

	// obtain RnaSeq metrics
	filtered_bams
		| get_metrics

	// Create scaled bigwig files for visualization
	CREATE_SCALED_BIGWIGS(filtered_bams)
}

workflow.onComplete {
	// Runs on success, cancel or fail
	WorkflowCompletion()
}

/*
 * fastqc
 * was originally part of extract_umi but the required container didn't support it
 */
process fastqc1 {
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path(fastqs)
	script:
		def r1 = fastqs.find { it.name.matches('.*[-_]R1[.-_].*') }
		def r2 = fastqs.find { it.name.matches('.*[-_]R2[.-_].*') }
		"""
          	#!/bin/bash

		module load fastqc

		echo "====== PROCESS_SUMMARY"
		echo "====== FASTQC1 ======"
		echo "Strategy: QC fastq files before extract_umi"
		echo " \$(fastqc --version)"
		echo "====== FASTQC1 ======"
		echo "====== PROCESS_SUMMARY"
		
		#QC
		fastqc -T ${task.cpus} -f fastq ${r1}
		fastqc -T ${task.cpus} -f fastq ${r2}
		"""
}

/*
 * fastqc
 */
process fastqc2 {
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path(fastqs)
	script:
		"""
          	#!/bin/bash

		module load fastqc

		echo "====== PROCESS_SUMMARY"
		echo "====== FASTQC2 ======"
		echo "Strategy: QC fastq files after extract_umi"
		echo " \$(fastqc --version)"
		echo "====== FASTQC2 ======"
		echo "====== PROCESS_SUMMARY"
		
		#QC
		fastqc -T ${task.cpus} -f fastq ${meta.id}_${meta.run}.extracted.R1.fastq.gz
		fastqc -T ${task.cpus} -f fastq ${meta.id}_${meta.run}.extracted.R2.fastq.gz
		"""
}


/*
 * extract umi
 * remove barcode and linker bases
 * append barcode to read name
 */
process extract_umi {
	container 'docker://quay.io/biocontainers/umi_tools:0.5.4--py27hdd9f355_1'
	params.extract_umi = [pattern: 'NNNNNNNNCCCCCC']
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path("${meta.id}_${meta.run}.extracted.*.fastq.gz")
	script:
		def r1 = fastqs.find { it.name.matches('.*[-_]R1[.-_].*') }
		def r2 = fastqs.find { it.name.matches('.*[-_]R2[.-_].*') }
		"""
          	#!/bin/bash

		echo "====== PROCESS_SUMMARY"
		echo "====== EXTRACT_UMI ======"
		echo "Strategy: Extract barcode and append to read name"
		echo " \$(umi_tools --version)"
		echo "pattern: ${params.extract_umi.pattern}"
		echo "N = UMI, C = barcode, X = reattached to read"
		echo "====== EXTRACT_UMI ======"
		echo "====== PROCESS_SUMMARY"
		
		# for kit SMART-Seq Total RNA Pico Input with UMIs (ZapR Mammalian)
		# UMI is the first 8 bases of read 2, followed by 6 linker bases to trim off
		# pattern = 'NNNNNNNNCCCCCC'

		echo "=== extracting barcode from ${meta.id}_${meta.run}"
				
		umi_tools extract \\
			-I ${r2} \\
			--read2-in=${r1} \\
			--bc-pattern=${params.extract_umi.pattern} \\
			--extract-method=string \\
			-S ${meta.id}_${meta.run}.extracted.R2.fastq.gz \\
			--read2-out=${meta.id}_${meta.run}.extracted.R1.fastq.gz \\
			--log=extract.log
		"""
}

/*
 * deduplicate with umi
 */
process dedup_umi {
	container 'docker://quay.io/biocontainers/umi_tools:0.5.4--py27hdd9f355_1'
	input:
		tuple val(meta), path(bam), path(bai)
	output:
		tuple val(meta), path("${meta.name}.umidedup.bam"), path(bam)
	script:
		"""
          	#!/bin/bash

		echo "====== PROCESS_SUMMARY"
		echo "====== DEDUP_UMI ======"
		echo "Strategy: Remove duplicates based on UMI"
		echo " \$(umi_tools --version)"
		echo "====== DEDUP_UMI ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== removing duplicated umi from ${meta.name}"

		umi_tools dedup -I $bam --paired --output-stats=deduplicated \\
			-S ${meta.name}.umidedup.bam > dedup.out 2>&1


		"""
}
/*
 * qc after umi
 * originally bundled with dedup_umi but modules and containers aren't compatible
 */
process dedup_qc {
	input:
		tuple val(meta), path(dedup_bam), path(bam)
	output:
		tuple val(meta), path(dedup_bam)
	script:
		"""
          	#!/bin/bash

		# load module
		module load samtools

		echo "====== PROCESS_SUMMARY"
		echo "====== DEDUP_QC ======"
		echo "Strategy: Count after DEDUP_UMI"
		echo "====== DEDUP_QC ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== deduplicated ${meta.name}"

		# read count summary
		initial_reads=\$(samtools view -c $bam)
		final_reads=\$(samtools view -c $dedup_bam)
		echo "Initial reads: \$initial_reads"
		echo "Reads retained: \$final_reads"
		echo "Total reads removed: \$((\$initial_reads - \$final_reads))"
		echo "Retention rate: \$(echo "scale=2; \$final_reads * 100 / \$initial_reads" | bc)%"

		"""
}


/*
 * alignment with star
 */
process star_align {
	params.star = [
		twopassMode: 'Basic',
		outSAMtype: 'BAM SortedByCoordinate',
		outBAMsortingBinsN: '100',
		quantMode: 'TranscriptomeSAM',
		outWigType: 'bedGraph',
		outWigStrand: 'Unstranded',
		outSAMunmapped: 'Within'

	]
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path("${meta.id}_${meta.run}.raw.bam")
	script:
		def r1 = fastqs.find { it.name.contains('.1.fq') }
		def r2 = fastqs.find { it.name.contains('.2.fq') }
		"""
		#!/bin/bash

		# load required modules
		module load star
		module load samtools
		module load rsem

		echo "====== PROCESS_SUMMARY"
		echo "====== STAR_ALIGN ======"
		echo "Strategy: Align paired-end reads with star"
		echo "\$(star --version | head -1 | awk '{print \$1,\$2,\$3}')"
		echo "Reference: ${params.transcriptome_index}"
		echo "RSEM Index: ${params.rsem_index}"
		echo "Parameters: ${params.star}"
		echo "====== STAR_ALIGN ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== aligning ${meta.id}_${meta.run}" 

		STAR --genomeDir ${params.transcriptome_index} \\
			--readFilesIn ${r1} ${r2} \\
			--runThreadN ${task.cpus} \\
			--twopassMode ${params.star.twopassMode} \\
			--outSAMtype ${params.star.outSAMtype} \\
			--outBAMsortingBinsN ${params.star.outBAMsortingBinsN} \\
			--quantMode ${params.star.quantMode} \\
			--outWigType ${params.star.outWigType} \\
			--outWigStrand ${params.star.outWigStrand} \\
			--outSAMunmapped ${params.star.outSAMunmapped} > star.out 2>&1

		tail -1 star.out

		mv Aligned.sortedByCoord.out.bam ${meta.id}_${meta.run}.raw.bam

		# rename for multiqc ID parsing
		mv Log.final.out ${meta.id}_${meta.run}.Log.final.out
	
		bam_count=\$(samtools view -c ${meta.id}_${meta.run}.raw.bam)
		echo "${meta.id}_${meta.run}.raw.bam: \$bam_count reads aligned"

		echo "=== estimating transcriptome alignment with rsem"
		rsem-calculate-expression --paired-end -p ${task.cpus} \\
			--alignments --strandedness reverse --no-bam-output \\
   			Aligned.toTranscriptome.out.bam ${params.rsem_index} ${meta.id}_${meta.run} > rsem.out 2>&1

		tail -5 rsem.out
		"""
}



/*
 * Count alignments overlapping feature (genes, peaks...)
 */
process count_features {
	publishDir 'output/counts', mode: 'copy'
	params.counts = '-s 2 --largestOverlap -p'
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path("${meta.name}.counts"), path("${meta.name}.biotypes")
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load subread  # for featureCounts
		
		echo "====== PROCESS_SUMMARY"
		echo "====== COUNT_FEATURES ======"
		echo "Strategy: Count reads overlapping features"
		echo "subread (featureCounts) \$(subread --version)"
		echo "Reference: ${params.gene_models}"
		echo "Parameters: ${params.counts}"
		echo "====== COUNT_FEATURES ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== Counting transcriped aligned reads for ${meta.name}"

		# Count reads overlapping peaks using featureCounts
				
		featureCounts -T ${task.cpus} ${params.counts} -a ${params.gene_models} \\
			-o ${meta.name}.counts $bam
		featureCounts -T ${task.cpus} ${params.counts} -a ${params.gene_models} \\
			-o ${meta.name}.biotypes -g gene_biotype $bam

		# Total counts from featureCounts method (sum column 7, skip header lines)
		featureCounts_total=\$(tail -n +3 ${meta.name}.counts | awk '{sum+=\$7} END {print sum}')
		echo "\$featureCounts_total total counts"
		"""
}

/*
 * Run a basic DESeq analysis
 * contains initial pipeline
 * can edit .Rmd in output folder and run with -resume to use edited parameters
 */
process run_DESeq {
	publishDir 'output', mode: 'copy'
	input:
		val(count_data)
	output:
		tuple path("*.html"), path("DESeq_results")
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load R
		
		echo "====== PROCESS_SUMMARY"
		echo "====== RUN_DESEQ ======"
		echo "Strategy: Generate basic DESeq2 analysis"
		echo "\$(R --version)"
		echo "edit analysis output/.Rmd and re-run with -resume"
		echo "====== RUN_DESEQ ======"
		echo "====== PROCESS_SUMMARY"

		if ! test -f ${launchDir}/output/DESeq_analysis.Rmd; then
  			cp ${projectDir}/bin/DESeq_analysis.Rmd ${launchDir}/output/DESeq_analysis.Rmd
		fi

		echo "=== Running DESeq analysis on collected count data"

		Rscript -e "rmarkdown::render('${launchDir}/output/DESeq_analysis.Rmd', \\
			params=list(collected_counts='${count_data}', database='${params.database}'))"
		"""
}

/*
 * Collect RNA-seq metrics using Picard
 */
process get_metrics {
	publishDir 'output/metrics', mode: 'copy'
	input:
		tuple val(meta), path(bam)
	output:
		path "${meta.name}.rna_metrics"
	script:
		"""
		#!/bin/bash

		# load required modules
		module load picard/2.22.0
		
		echo "====== PROCESS_SUMMARY"
		echo "====== GET_METRICS ======"
		echo "Strategy: Collect RNA-seq metrics using Picard CollectRnaSeqMetrics"
		echo "picard version 2.22.0"
		echo "Reference: ${params.database}"
		echo "====== GET_METRICS ======"
		echo "====== PROCESS_SUMMARY"
		
		echo "=== Collecting picard RNA-seq metrics for ${meta.name}"

		REF="\$(echo ${params.database}/*refflat)"
		INT="\$(echo ${params.database}/*interval)"
		
		/usr/bin/java -Xmx20G -jar \$PICARD CollectRnaSeqMetrics \\
			REF_FLAT=\$REF \\
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\
			RIBOSOMAL_INTERVALS=\$INT \\
			INPUT=$bam \\
			OUTPUT=${meta.name}.rna_metrics
		
		echo "Output: ${meta.name}.rna_metrics"
		"""
}


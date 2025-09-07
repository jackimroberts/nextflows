#!/usr/bin/env nextflow

// Include subworkflows and shared help
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { WorkflowCompletion } from '../subworkflows/workflow_complete.nf'
include { getSharedHelp } from '../modules/shared_help'

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
      		| filter { sample_id, sample_name, condition, extra_data, fastq_file ->
          		fastq_file.name.startsWith(sample_id + "_")
      		}
      		| map { sample_id, sample_name, condition, extra_data, fastq_file ->
          		def run_part = fastq_file.name.replaceAll("^${sample_id}_", "").replaceAll("[-_][LR][12][.-_].*", "")
          		def meta = [id: sample_id, name: sample_name, condition: condition, run: run_part, extra: extra_data]
          		[[meta.id, meta.run], meta, fastq_file]
      		}
      		| groupTuple()  // Group by the [id, run] tuple
      		| map { id_run_tuple, meta_list, fastq_files ->
          		[meta_list[0], fastq_files]
      		}
	// Trim and align individual sequencing runs
		| fastqc1
		| extract_umi
		| fastqc2
		| adapter_trim 
		| star_align
	// Combine samples by sample id
		| map { meta, bam_file ->
			[meta.id, meta, bam_file]
		}
		| groupTuple
		| map { sample_id, meta_list, bam_files ->
			[meta_list[0], bam_files]
	 	}
	// Separate samples with single or multiple sequencing runs
		| branch {
                        multiple: it[1].size() > 1
                        single: it[1].size() == 1
                }
                | set { bam_branches }

        // Only merge files from multiple runs
	// Avoids unnecessary SLURM jobs
        bam_branches.multiple 
		| merge_run_bams
        // Combine both streams back together
                | mix( bam_branches.single )
                | filter_bams
		| dedup_umi
		| set { filtered_bams }

	// Count features in each sample
	filtered_bams
		| count_features
		| collect
		| run_DESeq

	// obtain RnaSeq metrics
	filtered_bams
		| get_metrics

	// Calculate read depth for each bam
	filtered_bams
		| measure_depth
		| map { meta, bam_file, depth ->
          		def depth_value = depth.trim().split('\n').last() as Integer
   
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
		fastqc -T \$(nproc) -f fastq ${r1}
		#fastqc -T \$(nproc) -f fastq ${r2}
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
		#fastqc -T \$(nproc) -f fastq ${meta.id}_${meta.run}.extracted.R1.fastq.gz
		#fastqc -T \$(nproc) -f fastq ${meta.id}_${meta.run}.extracted.R2.fastq.gz
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
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path("${meta.name}.umidedup.bam")
	script:
		"""
          	#!/bin/bash

		# load module
		module load umi-tools
		#module load samtools

		echo "====== PROCESS_SUMMARY"
		echo "====== DEDUP_UMI ======"
		echo "Strategy: Remove duplicates based on UMI"
		echo " \$(umi_tools --version)"
		echo "====== DEDUP_UMI ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== removing duplicated umi from ${meta.name}"

		umi_tools dedup -stdin=$bam --paired --output-stats=deduplicated \\
			-S ${meta.name}.umidedup.bam

		# read count summary
		#initial_reads=\$(samtools view -c $bam)
		#final_reads=\$(samtools view -c ${meta.name}.umidedup.bam)
		#echo "Initial reads: \$initial_reads"
		#echo "Reads retained: \$final_reads"
		#echo "Total reads removed: \$((\$initial_reads - \$final_reads))"
		#echo "Retention rate: \$(echo "scale=2; \$final_reads * 100 / \$initial_reads" | bc)%"

		"""
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

		echo "====== PROCESS_SUMMARY"
		echo "====== ADAPTER_TRIM ======"
		echo "Strategy: Remove adapters and low quality bases from reads"
		echo "trimmed with cutadapt \$(cutadapt --version)"
		echo "-O 1 --nextseq-trim=20 -m 1"
		echo "forward adapter: \$adaptf"
		echo "and reverse: \$adaptr"
		echo "====== ADAPTER_TRIM ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== trimming ${meta.id}_${meta.run}"

		cutadapt -O 1 --nextseq-trim=20 -m 1 -a \$adaptf -A \$adaptr -j \$(nproc) \\
			-o ${meta.id}_${meta.run}.1.fq -p ${meta.id}_${meta.run}.2.fq \\
			$fastqs > cutadapt.log
		
		# Extract key stats from cutadapt log
		sed -n '/Total read pairs processed/,/=== First read: Adapter 1 ===/p' cutadapt.log | head -n -1
				
		#QC
		fastqc -T \$(`nproc`) -f fastq ${meta.id}_${meta.run}.1.fq
		fastqc -T \$(`nproc`) -f fastq ${meta.id}_${meta.run}.2.fq

		# %rRNA and species using fastq screen
		${params.fastq_screen} --conf ${params.fastq_screen_config} \\
			--threads \$(nproc) --subset 1000000 ${meta.id}_${meta.run}.1.fq
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
		echo "Parameters: ${params.star}"
		echo "====== STAR_ALIGN ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== aligning ${meta.id}_${meta.run}" 

		STAR --genomeDir ${params.transcriptome_index} \\
			--readFilesIn ${r1} ${r2} \\
			--runThreadN \$(`nproc`) \\
			--twopassMode ${params.star.twopassMode} \\
			--outSAMtype ${params.star.outSAMtype} \\
			--outBAMsortingBinsN ${params.star.outBAMsortingBinsN} \\
			--quantMode ${params.star.quantMode} \\
			--outWigType ${params.star.outWigType} \\
			--outWigStrand ${params.star.outWigStrand} \\
			--outSAMunmapped ${params.star.outSAMunmapped} \\
			> star.log 2>&1
		
		# Output only the last line of STAR stdout
		tail -1 star.log

		mv Aligned.sortedByCoord.out.bam ${meta.id}_${meta.run}.raw.bam

		# rename for multiqc ID parsing
		mv Log.final.out ${meta.id}_${meta.run}.Log.final.out
	
		bam_count=\$(samtools view -c ${meta.id}_${meta.run}.raw.bam)
		echo "${meta.id}_${meta.run}.raw.bam: \$bam_count reads aligned"

		# QC
		rsem-calculate-expression --paired-end -p \$(`nproc`) \\
			--alignments --strandedness reverse --no-bam-output \\
   			Aligned.toTranscriptome.out.bam ${params.rsem_index} ${meta.id}_${meta.run}
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
		
		echo "====== PROCESS_SUMMARY"
		echo "====== MERGE_RUN_BAMS ======"
		echo "Strategy: Combine BAM files from the same sample across runs"
		echo "\$(samtools --version | head -1)"
		echo "====== MERGE_RUN_BAMS ======"
		echo "====== PROCESS_SUMMARY"
		

		echo "=== Merging BAM files from multiple runs"
		for bam in $bam_files; do
			read_count=\$(samtools view -c \$bam)
			echo "\$bam: \$read_count aligned reads"
		done

		samtools merge ${meta.id}.merged.bam $bam_files

		echo "=== Merged bam"
		merged_count=\$(samtools view -c ${meta.id}.merged.bam)
		echo "${meta.id}.merged.bam: \$merged_count aligned reads"

		"""
}

/*
 * filter bams for duplicates, unaligned, unmated, mitochondria, blacklist
 */
process filter_bams {
	params.filter = [flag: 3340]
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path('*.filtered.bam')
	script:
		"""
		sh ${projectDir}/bin/filter_bams.sh $bam ${meta.name} \\
			${params.genome_blacklist} ${params.filter.flag}
		"""
}



/*
 * Count alignments overlapping feature (genes, peaks...)
 */
process count_features {
	publishDir 'output/counts', mode: 'copy'
	params.counts = '--stranded reverse --largestOverlap --paired'
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
				
		featureCounts -T \$(`nproc`) ${params.counts} -a ${params.gene_models} \\
			-o ${meta.name}.counts $bam
		featureCounts -T \$(`nproc`) ${params.counts} -a ${params.gene_models} \\
			-o ${meta.name}.biotypes -g gene_biotype $bam

		# Total counts from featureCounts method (sum column 7, skip header lines)
		featureCounts_total=\$(tail -n +3 ${meta.name}.counts | awk '{sum+=\$7} END {print sum}')
		echo "\$featureCounts_total total counts"
		"""
}

/*
 * Run a basic DESeq analysis
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
		echo "====== RUN_DESEQ ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== Running DESeq analysis on collected count data"

		Rscript -e "rmarkdown::render('${projectDir}/bin/DESeq_analysis.Rmd', \\
			params=list(count_data='${count_data}', database='${params.database}'))"
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
		
		echo "====== PROCESS_SUMMARY"
		echo "====== MEASURE_DEPTH ======"
		echo "Strategy: Count read pairs for scaling normalization"
		echo "\$(samtools --version | head -1)"
		echo "====== MEASURE_DEPTH ======"
		echo "====== PROCESS_SUMMARY"
		
		echo "$bam depth:"
		## Calculate size for eventual scaling
		samtools view -f 64 -c $bam

		"""
}

/*
 * convert bam files to bigwig
 */
process bam_to_bigwig {
	publishDir 'output/bigwigs', mode: 'copy'
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

		echo "====== PROCESS_SUMMARY"
		echo "====== BAM_TO_BIGWIG ======"
		echo "Strategy: Convert BAM to normalized BigWig for visualization"
		echo "deeptools \$(deeptools --version)"
		echo "\$(samtools --version | head -1)"
		echo "Scale factor: min_depth/sample_depth"
		echo "====== BAM_TO_BIGWIG ======"
		echo "====== PROCESS_SUMMARY"

		# create scale based on depth and minimum depth of all samples
		scale=\$(echo "scale=5; $min_depth/${meta.depth}" | bc)

		# index bam and generate bigwig
		samtools index $bam
		bamCoverage -b $bam --scaleFactor \$scale -o ${meta.name}.bw
		
		file_size=\$(du -h "${meta.name}.bw" | cut -f1)
		echo "Created ${meta.name}.bw: \$file_size"
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
		
		echo "=== Collecting RNA-seq metrics for ${meta.name}"

		REF="\$(echo ${params.database}/*refflat)"
		INT="\$(echo ${params.database}/*interval)"
		
		/usr/bin/java -Xmx20G -jar \$PICARD CollectRnaSeqMetrics \\
			REF_FLAT=\$REF \\
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\
			RIBOSOMAL_INTERVALS=\$INT \\
			INPUT=$bam \\
			OUTPUT=${meta.name}.rna_metrics
		
		echo "Metrics Collection Complete"
		echo "Output: ${meta.name}.rna_metrics"
		"""
}


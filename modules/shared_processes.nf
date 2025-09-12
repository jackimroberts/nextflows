/*
 * Shared processes for ATAC-seq and scRNA-seq pipelines
 */

/*
 * generate table with specific format:
 * sampleID sampleName expCondition
 */
process make_sample_table {
	input:
		path input_sample_table
	output:
		path "sample_table.tsv"
	script:
		"""
		module load R

		echo "====== PROCESS_SUMMARY"
		echo "====== MAKE_SAMPLE_TABLE ======"
		echo "Strategy: Converts txt to .tsv [id, name, condition ...]"
		echo "Manipulation in \$(R --version | head -1)"
		echo "====== MAKE_SAMPLE_TABLE ======"
		echo "====== PROCESS_SUMMARY"

		Rscript ${projectDir}/../bin/make_sample_table.R $input_sample_table
		"""
}

/*
 * Pull fastq files
 */
process get_fastq {
	publishDir 'output', pattern: '**md5*', mode: 'copy', overwrite: true
	input:
		val input_file_source
		path sample_table_file
	output:
		path "**{fastq,md5}*"
	script:
		"""
		echo "====== PROCESS_SUMMARY"
		echo "====== GET_FASTQ ======"
		echo "Strategy: Downloads fastq files from storage"
		echo "====== GET_FASTQ ======"
		echo "====== PROCESS_SUMMARY"

		sh ${projectDir}/../bin/get_files.sh "${input_file_source}" "${sample_table_file}" "${launchDir}" "${task.cpus}"
		"""
}

/*
 * decompress ora files to *gz format
 */
process decompress {
	input:
		path ora_file
	output:
		path "*fastq.gz"
	script:
		"""
		# module load oradecompression/2.7.0

		echo "====== PROCESS_SUMMARY"
		echo "====== DECOMPRESS ======"
		echo "Strategy: Decompress .ora files to .fastq.gz"
		echo "Uses orad 2.7.0"
		echo "====== DECOMPRESS ======"
		echo "====== PROCESS_SUMMARY"
		
		echo "=== input ORA file"
		ls -lh ${ora_file} | awk '{print \$9 ": " \$5}'

		/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/orad/2.7.0/orad "${ora_file}" >/dev/null 2>&1

		echo "=== decompressed fastq.gz files"
		ls -lh *.fastq *.fastq.gz | awk '{print \$9 ": " \$5}'

		"""
}

/*
 * Make mini fastq for testing
 */
process miniaturize {
	stageInMode 'copy'
	input:
		path fastq_file
	output:
		path "${fastq_file}"
	script:
		"""
		echo "====== PROCESS_SUMMARY"
		echo "====== MINIATURIZE ======"
		echo "Strategy: Miniaturize fastq files for quicker pipelines"
		echo "Files reduced and compression is maintained"
		echo "====== MINIATURIZE ======"
		echo "====== PROCESS_SUMMARY"

		echo "=== input file"
		ls -lh *.fastq *.fastq.gz | awk '{print \$9 ": " \$5}'

        	if [[ "${fastq_file}" == *.gz ]]; then
            		zcat $fastq_file | head -n 1000000 | gzip > temp_${fastq_file}
        	else
            		head -n 1000000 $fastq_file > temp_${fastq_file}
        	fi
		mv temp_${fastq_file} ${fastq_file}

		echo "=== output file"
		ls -lh *.fastq *.fastq.gz | awk '{print \$9 ": " \$5}'

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

		cutadapt -O 1 --nextseq-trim=20 -m 1 -a \$adaptf -A \$adaptr -j ${task.cpus} \\
			-o ${meta.id}_${meta.run}.1.fq -p ${meta.id}_${meta.run}.2.fq \\
			$fastqs > ${meta.id}_${meta.run}_cutadapt.log
		
		# Extract key stats from cutadapt log
		sed -n '/Total read pairs processed/,/=== First read: Adapter 1 ===/p' ${meta.id}_${meta.run}_cutadapt.log | head -n -1
		
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
 * filter bams for duplicates, unaligned, mitochondria, blacklist
 */
process filter_bams {
	input:
		tuple val(meta), path(bam)
	output:
		tuple val(meta), path('*.filtered.bam'), path('*.filtered.bam.bai')
	script:
		"""
		sh ${projectDir}/../bin/filter_bams.sh $bam ${meta.name} \\
			${params.genome_blacklist} ${params.filter.flag} ${task.cpus}
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
 * Collect QC files from successful processes
 */
process collect_qc {
	publishDir "${params.outputDir}/qc_metrics", mode: 'move'
	input:
		val trigger
	output:
		path "qc_metrics"
	script:
		"""
		#!/bin/bash
		
		echo "====== PROCESS_SUMMARY"
		echo "====== COLLECT_QC ======"
		echo "Strategy: Collect QC files from completed processes"
		echo "Searches work directories for QC log files"
		echo "====== COLLECT_QC ======"
		echo "====== PROCESS_SUMMARY"
		
		mkdir -p qc_metrics
		
		# Collect cutadapt logs
		find ${launchDir}/work/*/*/*cutadapt.log -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect STAR logs
		find ${launchDir}/work/*/*/*Log.final.out -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect MACS2 logs
		find ${launchDir}/work/*/*/*macs.log -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect Bowtie2 logs
		find ${launchDir}/work/*/*/*bowtie.log -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect FastQC files
		find ${launchDir}/work/*/*/*_fastqc.html -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		find ${launchDir}/work/*/*/*_fastqc.zip -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect featureCounts summaries
		find ${launchDir}/work/*/*/*summary -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect STAR stat files
		find ${launchDir}/work/*/*/*stat -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Collect RNA metrics files
		find ${launchDir}/work/*/*/*rna_metrics -exec cp {} qc_metrics/ \\; 2>/dev/null || true
		
		# Show what was collected
		echo "=== Collected QC files:"
		ls -la qc_metrics/ || echo "No QC files found"
		"""
}

/*
 * MultiQC - aggregate all QC reports
 */
process multiqc {
	publishDir "${params.outputDir}", mode: 'copy'
	input:
		path(qc_metrics_dir)
	output:
		path "multiqc_report.html"
		path "multiqc_data"
	script:
		"""
		#!/bin/bash
		
		# Load required modules
		module load multiqc
		
		echo "====== PROCESS_SUMMARY"
		echo "====== MULTIQC ======"
		echo "Strategy: Aggregate all QC reports into single report"
		echo "\$(multiqc --version)"
		echo "====== MULTIQC ======"
		echo "====== PROCESS_SUMMARY"
		
		echo "=== Generating MultiQC report from ${qc_metrics_dir}"
		
		multiqc ${qc_metrics_dir} --filename multiqc_report.html
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
		echo "Created ${meta.name}.bw, file size: \$file_size"
		"""
}
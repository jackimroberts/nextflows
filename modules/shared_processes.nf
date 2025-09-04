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

		echo "=== PROCESS_MAKE_SAMPLE_TABLE ==="
		echo "Strategy: Converts *.txt to tsv [id, name, condition]"
		echo "Manipulation in R \$(R --version)"
		echo "=== PROCESS_MAKE_SAMPLE_TABLE ==="

		Rscript ${projectDir}/../shared_bin/make_sample_table.R $input_sample_table
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
		echo "=== PROCESS_GET_FASTQ ==="
		echo "Strategy: Downloads fastq files from storage"
		echo "=== PROCESS_GET_FASTQ ==="

		sh ${projectDir}/../shared_bin/get_files.sh "${input_file_source}" "${sample_table_file}" "${launchDir}"
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
		echo "=== PROCESS_DECOMPRESS ==="
		echo "Strategy: Decompress *ora files to *fastq.gz"
		echo "Uses orad 2.7.0"
		echo "=== PROCESS_DECOMPRESS ==="
		
		echo "=== input ORA file"
		du -h ${ora_file}
		
		# module load oradecompression/2.7.0
		/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/orad/2.7.0/orad "${ora_file}"

		echo "=== decompressed fastq.gz files"
		ls -lh *.fastq *.fastq.gz | awk '{print $9 ": " $5}'

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
		echo "=== PROCESS_MINIATURIZE ==="
		echo "Strategy: Miniaturize fastq files for quicker pipelines
		echo "If *gz it decompresses and recompresses
		echo "=== PROCESS_MINATURIZE ==="

		echo "=== input file"
		ls -lh *.fastq *.fastq.gz | awk '{print $9 ": " $5}'

        	if [[ "${fastq_file}" == *.gz ]]; then
            		zcat $fastq_file | head -n 1000000 | gzip > temp_${fastq_file}
        	else
            		head -n 1000000 $fastq_file > temp_${fastq_file}
        	fi
		mv temp_${fastq_file} ${fastq_file}

		echo "=== output file"
		ls -lh *.fastq *.fastq.gz | awk '{print $9 ": " $5}'

		"""
}
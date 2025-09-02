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
		sh ${projectDir}/../shared_bin/get_files.sh "${input_file_source}" "${sample_table_file}" "${launchDir}"
		"""
}

/*
 * decompress ora files to *gz format
 */
process decompress {
	input:
		path fastq_file
	output:
		path "*fastq.gz"
	script:
		"""
		if [[ "${fastq_file}" == *".ora" ]]; then
			# module load oradecompression/2.7.0
			/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/orad/2.7.0/orad "${fastq_file}" --rm
		fi
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
        	if [[ "${fastq_file}" == *.gz ]]; then
            		zcat $fastq_file | head -n 1000000 | gzip > temp_${fastq_file}
        	else
            		head -n 1000000 $fastq_file > temp_${fastq_file}
        	fi
		mv temp_${fastq_file} ${fastq_file}
		"""
}
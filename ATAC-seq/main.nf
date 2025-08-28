#!/usr/bin/env nextflow

params.miniaturize = false
params.fastq_source = true
params.sample_table = "*.txt"

/*
 * Things I'd like to add:
 ** Ability to get SRA files
 ** add sample names and conditions
 ** tidy up all scripts with headers
 ** qc throughout
 ** create container
 ** deploy on cloud
*/

params.help = false
if (params.help) {
        println """
        ATAC-seq nextflow pipeline

	--sample_table *txt
		Must be in launch folder. Table from gnomex

	--fastq_source="java -jar ./fdt....gnomex..."
			Pulls fastq files from gnomex
				get this command from:
				gnomex > navigate to experiment > "Files" > "Download Files" > 
					move fastq folder to the right > "FDT Command Line" > copy command
		="SSD/YYYYMMDD_run_identifier/email_subject_line:password"
			UCSF core emails a filepath and password after sequencing 
		="SRA"
			Using sample table text file, downloads fastqs from SRA repository
			One SRA ID per line (SRR...), skips comments and empty lines
			(default: true)
	
	--miniaturize true/false
		Create mini fastq files of 2,500 reads for testing (default: false)


        """
        exit 0
    }

workflow {

	// wrangle sample sheet into specific format: id, name, condition
	Channel.fromPath(params.sample_table)
		| make_sample_table 
		| splitCsv(header:['sample_id','sample_name','condition'], sep:"\t") 
		| set {sample_sheet}
      	
	// Find existing fastq files
	Channel.fromPath("**fastq*", checkIfExists: false)
		| filter { it.exists() && !it.toString().contains('/work/') } // exclude those in 'work'
		| set { existing_files }

	// Downloaded files if specified and combine with existing files
	// Separate files by type and handle appropriately
	(params.fastq_source != true ? Channel.value(params.fastq_source) | get_fastq : Channel.empty())
		| mix(existing_files)
		| branch {
			ora_files: it.name.endsWith('.ora')
			gz_files: it.name.endsWith('.gz')
		}

	// Only decompress .ora files, then mix with .gz files
	ora_files
		| decompress
		| mix(gz_files)
		| flatten
		| set { flat_fastq_list }

	// Downsize fastq for faster runs
        if (params.miniaturize == true) {
            flat_fastq_list
		| miniaturize
		| set { fastq_list }
        } else {
            flat_fastq_list
		| set { fastq_list }
        }
		
	// pair up fastq files by sample_id
	fastq_list
		| map{ file ->
            		def sample_id = file.name.replaceAll("_.*", "")
            		tuple(sample_id, file)
        	  }
        	| groupTuple
		| set {paired_fastq}
	
	// run pipeline for individual samples
	paired_fastq 
		| adapter_trim 
		| bowtie_align 
		| filter_bams
	
	// replace sample number with name
	// keep name, bamfile, depth for bam to bw
	filter_bams.out
		| join (sample_sheet)
		| map { row -> [row[3], row[1], row[2]] }
		| set {named_filtered_bams}
	view(named_filtered_bams)
	// find min_depth and append to tuple
	named_filtered_bams
		| map{ row -> row[2]}
		| min()
		| combine(named_filtered_bams) 
		| bam_to_bigwig
	
	// replace sample number with name
	// keep name, bamfile, condition for peakcalling
	filter_bams.out
		| join (sample_sheet)
		| map { row -> [row[3], row[1], row[4]] }
		| set {bam_sample_sheet} 
	
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
/*
	// run Tim's pipeline
	bam_sample_sheet
		| map{ row -> [row[2],row[1]] }
		| groupTuple
		| tuple_to_stdout  
		| collect
		| create_multimacs_run 
		//| run_multimacs | view
*/
}

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
		Rscript ${projectDir}/bin/make_sample_table.R $input_sample_table
		"""
}

/*
 * Pull fastq files
 */
process get_fastq {
	input:
		val input_file_source
	output:
		path "**fastq*"
	script:
		"""
		sh ${projectDir}/bin/get_files.sh "${input_file_source}"
		"""
}

/*
 * decompress ora files to *gz format
 */
process decompress {
	input:
		path fastq_file
	output:
		path "*fastq.gz", includeInputs: true
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
		head -n 100000 $fastq_file > temp_${fastq_file}
		mv temp_${fastq_file} ${fastq_file}
		"""
}

/*
 * adapter trimming
 */
process adapter_trim {
	input:
		tuple val(sample_num), path(fastqs)
	output:
		tuple val(sample_num), path('*fq')
	script:
		"""
		sh ${projectDir}/bin/adapter_trim.sh $sample_num $fastqs
		"""
}

/*
 * alignment with bowtie2, using mm10
 */
process bowtie_align {
	input:
		tuple val(sample_num), path(fastqs)
	output:
		tuple val(sample_num), path('*raw.bam')
	script:
		"""
		sh ${projectDir}/bin/bowtie_align.sh $fastqs $sample_num ${params.genome_index}
		"""
}

/*
 * filter bams for duplicates, unaligned, mitochondria, blacklist
 * also append sequencing depth to output
 */
process filter_bams {
	input:
		tuple val(sample_num), path(bam)
	output:
		tuple val(sample_num), path('*filtered.bam'), stdout
	script:
		"""
		sh ${projectDir}/bin/filter_bams.sh $bam $sample_num ${params.genome_blacklist}
		"""
}

/*
 * convert bam files to bigwig
 */
process bam_to_bigwig {
	publishDir 'bigwigs', mode: 'copy'
	input:
		tuple val(min_depth), val(sample_name), path(bam), val(depth)
	output:
		path "*bw"
	script:
		"""
		sh ${projectDir}/bin/bam_to_bigwig.sh $sample_name $bam $depth $min_depth
		"""
}

/*
 * run filtered bam files through Tim's pipeline
 */
process tuple_to_stdout {
	input: 
		tuple val(condition), val(bam_files)
	output: 
		stdout
	script:
		"""
		echo --chip $bam_files
		echo --name $condition
		"""
}
process create_multimacs_run {
publishDir 'multimacs_pipeline', mode: 'copy'
	input:
		stdin str
	output:
		path "*sh"
	script:
		"""
		sh ${projectDir}/bin/create_multimacs_run.sh ${projectDir}/bin/multimacs_run_template.txt
		"""
}
process run_multimacs {
	publishDir 'multimacs_pipeline', mode: 'copy'
	input:
		path multimacs_command
	output:
		path "*"
	script:
		"""
		sh $multimacs_command
		"""
}



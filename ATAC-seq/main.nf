#!/usr/bin/env nextflow

params.mini = true
params.fastq_source = "empty"
params.sample_table = "*.txt"

/*
 * Things I'd like to add:
 ** Ability to get SRA files
 ** change publishDir based on sample names
 ** add sample names and conditions
 ** tidy up all scripts with headers
 ** qc throughout
 ** upload to github
 ** create container
 ** deploy on cloud
*/

params.help = false
if (params.help) {
        println """
        This pipeline is designed to process ATAC-seq data beginning from fastq files

	--sample_table *txt
		Must be in this folder. Table from gnomex

	--fastq_source "java -jar ./fdt....gnomex..."
		Currently broken
		Pulls fastq files from gnomex
		get this command from:
		gnomex > navigate to experiment > "Files" > "Download Files" > 
			move over fastq folder > "FDT Command Line" > copy command
	
	--mini	subsets fastq files to 100k reads

        """
        exit 0
    }

input_fastq_source = Channel.value(params.fastq_source)
input_sample_details = Channel.fromPath(params.sample_table)

workflow {

	// wrangle sample sheet into specific format: id, name, condition
	input_sample_details 
		| make_sample_table 
		| splitCsv(header:['sample_id','sample_name','condition'], sep:"\t") 
		| set {sample_sheet}

	// get fastq files organized
	//input_fastq_source 
	//	| get_fastq
	// broken. will run concurrently with next step. Requires piping together, conditonal statement or completed signal

	Channel.fromPath("**fastq.gz")
		| filter { !it.toString().contains('/work/') }
		| flatten 
		| set {fastq_list}
	
	// make fastq and sample sheet 'mini' if specified
 	if (params.mini) { // helpful for testing the pipeline quickly
		// adjust sample names for clarity
		sample_sheet
			| map{ row -> 
				["mini.${row.sample_id}","mini.${row.sample_name}",row.condition] 
			  }
			| set {sample_sheet}
		// subset fastqs and rename
		fastq_list 
			| make_mini 
			| set {fastq_list}
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
	// run Tim's pipeline
	bam_sample_sheet
		| map{ row -> [row[2],row[1]] }
		| groupTuple
		| tuple_to_stdout  
		| collect
		| create_multimacs_run 
		//| run_multimacs | view
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
		path "**fastq.gz"
	script:
		"""
		sh ${projectDir}/bin/get_files.sh "${input_file_source}"
		"""
}

/*
 * Make mini fastq for testing
 */
process make_mini {
	input:
		path fastq_files
	output:
		path "mini*"
	script:
		"""
		head -n 100000 $fastq_files > mini.$fastq_files
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



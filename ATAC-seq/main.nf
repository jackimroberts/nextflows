#!/usr/bin/env nextflow

params.miniaturize = false
params.fastq_source = true
params.sample_table = "*.txt"

/*
 * Things I'd like to add:
 ** tidy up all scripts with headers
 ** qc throughout
 ** create containerized version
 ** deploy on cloud
*/

params.help = false
if (params.help) {
        println """
        ATAC-seq nextflow pipeline

        REQUIRED PARAMETERS:
        --sample_table *txt
                Must be in launch folder. Table from gnomex or tsv format: id name condition

        OPTIONAL PARAMETERS:
        --miniaturize true/false
                Create mini fastq files of 2,500 reads for testing (default: false)

        --fastq_source="source1,source2,..."
                Specify input sources, can combine multiple with commas
                
                ="java -jar ./fdt....gnomex..."
                        Pulls fastq files from gnomex
                        Get this command from:
                        gnomex > navigate to experiment > "Files" > "Download Files" > 
                        move fastq folder to the right > "FDT Command Line" > copy command
                
                ="CoreBrowser"
                        Retrieves unarchived files from Utah Core Browser
                        Select files > More options dropdown menu > Secure link >
                        paste into a "core_links" text file
                
                ="SSD/YYYYMMDD_run_identifier/email_subject_line:password"
                        UCSF core emails a filepath and password after sequencing
                
                ="SRA"
                        Using sample table text file, downloads fastqs from SRA repository
                        Finds any SRR#+, skipping comment and empty lines
                
                (default: searches for existing fastq files in launch directory)

        """
        exit 0
}

workflow {

	// wrangle sample sheet into specific format: id, name, condition
	Channel.fromPath(params.sample_table)
		| make_sample_table 
		| splitCsv(header:['sample_id','sample_name','condition'], sep:"\t") 
		| map { row -> [row.sample_id, row.sample_name, row.condition] }
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
			ora: it.name.endsWith('.ora')
			other: true // should catch *fastq and *fastq.gz
		}
		| set { files }

	// Only decompress .ora files, then mix with .gz files
	files.ora
		| decompress
		| mix(files.other)
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
			| map { meta, bam_file ->
				[meta.id, meta, bam_file]
		}
		| groupTuple
	// Combine samples from different runs
	// Merge bam files, filter, and calculate depth
		| map { sample_id, meta_list, bam_files ->
			[meta_list[0], bam_files]
	 	}
		| merge_run_bams
		| filter_bams
		| map { meta, bam_file, depth ->
          		def more_meta = meta + [depth: depth.trim()]
          		[more_meta, bam_file]
		}
		| set { filtered_bams_with_depth }

  	// Find minimum depth across all samples for scaling
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
	publishDir 'output', pattern: '**md5*', mode: 'copy', overwrite: true
	input:
		val input_file_source
	output:
		path "**fastq*"
	script:
		"""
		sh ${projectDir}/bin/get_files.sh "${input_file_source}" "${params.sample_table}"
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
		head -n 100000 $fastq_file > temp_${fastq_file}
		mv temp_${fastq_file} ${fastq_file}
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
		sh ${projectDir}/bin/adapter_trim.sh ${meta.id}_${meta.run} $fastqs
		"""
}

/*
 * alignment with bowtie2, using mm10
 */
process bowtie_align {
	input:
		tuple val(meta), path(fastqs)
	output:
		tuple val(meta), path('*raw.bam')
	script:
		"""
		sh ${projectDir}/bin/bowtie_align.sh $fastqs ${meta.id}_${meta.run} ${params.genome_index}
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
            ln -s $bam_files ${meta.id}.merged.bam
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
		sh ${projectDir}/bin/bam_to_bigwig.sh ${meta.name} $bam ${meta.depth} $min_depth
		"""
}

/*
 * run filtered bam files through Tim's pipeline
 */
process tuple_to_stdout {
	input: 
		tuple val(meta), val(bam_files)
	output: 
		stdout
	script:
		"""
		echo --chip $bam_files
		echo --name ${meta.condition}
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



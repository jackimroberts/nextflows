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
		| flatten
		| branch {
			ora: it.name.endsWith('.ora')
			gz: it.name.endsWith('.gz') || it.name.endsWith('.fastq')
		}
		| set { files }

	// Only decompress .ora files, then mix with .gz files
	files.ora
		| decompress
		| mix(files.gz)
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
		//| collect
		//| run_multimacs
}

workflow.onComplete {
    // Run resource usage analysis
    def workDir = workflow.launchDir
    def analyzerScript = "${projectDir}/../shared_bin/slurm_usage_analyzer.sh"
    
    println """
    ==========================================================
    Pipeline execution summary
    ==========================================================
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work directory: ${workDir}
    ==========================================================
    """
    
    if (file(analyzerScript).exists()) {
        println "Running resource usage analysis..."
        try {
            def proc = ["bash", analyzerScript, workDir].execute()
            proc.waitFor()
            
            if (proc.exitValue() == 0) {
                println proc.text
            } else {
                println "Resource analysis completed with warnings:"
                println proc.text
                if (proc.err.text) {
                    println "Error output:"
                    println proc.err.text
                }
            }
        } catch (Exception e) {
            println "Failed to run resource analysis: ${e.message}"
        }
    } else {
        println "Resource analyzer script not found at: ${analyzerScript}"
        println "Skipping resource usage analysis."
    }
    
    println "=========================================================="
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
	output:
		path "**fastq*"
	script:
		"""
		sh ${projectDir}/../shared_bin/get_files.sh "${input_file_source}" "${params.sample_table}" "${launchDir}"
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



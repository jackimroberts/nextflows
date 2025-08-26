#!/usr/bin/env nextflow

params.expected_cell_number = 10000
params.fastq_source = true
//params.fastq_source = "empty"
params.sample_table = "*txt"

params.help = false
if (params.help) {
        println """
        This pipeline is designed to pull fastq files from Gnomex and run scRNA-seq

	--sample_table *txt
		Must be in this folder. 
		Table from gnomex contains "ID" and "Sample Name"
			gnomex > navigate to experiment > "Experiment Design" > "Download Sample Sheet" > 
				move to the launch folder
		Can be manually generated with format "ID" "name" "condition"
			ID is the unique fastq file prefix

	--fastq_source "java -jar ./fdt....gnomex..."
		Pulls fastq files from gnomex
			get this command from:
			gnomex > navigate to experiment > "Files" > "Download Files" > 
				move fastq folder to the right > "FDT Command Line" > copy command
		default is 'true' and will search for **fastq.gz in launch folder
	
	--expected_cell_number 10000
		Used by CellRanger default 10,000 cells	

        """
        exit 0
    }

workflow {

	// wrangle sample sheet into specific format: id, name, condition
	Channel.fromPath(params.sample_table)
		| make_sample_table
		| splitCsv( sep:"\t") 
		| set {sample_sheet}

	// Download fastq files
	//input_fastq_source 
Channel.value(params.fastq_source)
		| get_fastq // download if command is given
      	
	// Find existing fastq files
	Channel.fromPath("**/*fastq*", checkIfExists: false)
		| filter { it.exists() && !it.toString().contains('/work/') } // exclude those in 'work'
		| set { existing_files }

	// Combine downloaded and existing files
	get_fastq.out
		.mix(existing_files)
		| decompress // Handles *ora compression
		| flatten // flattens list
		| set set { fastq_list }

	// pair up fastq files by sample id
	// replace sample id with name
	fastq_list
		| map{ file ->
            		def sample_id = file.name.replaceAll("_.*", "")
            		tuple(sample_id, file)
        		}
        		| groupTuple
		| join (sample_sheet)
		| map { row -> [row[0], row[2], row[1]] }
		| set {paired_fastq}

	// run pipeline for individual samples
	paired_fastq 
		| run_cellranger
		| collect 
		| seurat_markdown
		//| make_shiny

/*
	// run pipeline for individual samples
	paired_fastq 
		| run_cellranger
		| run_velocyto

	run_velocyto.out.file_tree 
		| collect 
		| seurat_markdown
		//| make_shiny
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
		path "**/*fastq*"
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
process make_mini {
	input:
		path fastq_files
	output:
		path "mini*"
	script:
		"""
		head -n 1000000 $fastq_files > mini$fastq_files
		"""
}

/*
 * run cellRanger
 */
process run_cellranger {
	//publishDir 'output', pattern: '**html', mode: 'copy', overwrite: true
	input:
		tuple val(sample_id), val(sample_name), path(fastqs)
	output:
		path "${sample_name}"
	script:
		"""
		#place sample fastqs in their individual folder
		mkdir ${sample_id}_fastqs
		cp $fastqs ${sample_id}_fastqs
		
		#submit to cellranger
		sh ${projectDir}/bin/run_cellranger.sh \
			$sample_id \
			${sample_id}_fastqs \
			${params.cellRanger_transcriptome} \
			${params.expected_cell_number}

		#change id to name from here
		mv ${sample_id} ${sample_name}
		rm -r ${sample_id}_fastqs
		"""
}
/*
* rm -r ${sample_id}_fastqs
*/
		
/*
 * run velocyto
 * produces loom file containing spliced and unspliced information
 */
process run_velocyto {
	publishDir 'output', pattern: '**html', mode: 'copy', overwrite: true
	//runs out of memory in parallel
	maxForks 1
	input:
		path cellRanger_out
	output:
		path "${cellRanger_out}", emit: file_tree
		path '**html', emit: loom_file
	script:
		"""
		module load velocyto
		module load samtools
		
		which velocyto
		which samtools
	
		if ! test -f ${launchDir}/work/genes.gtf; then
			cp ${params.cellRanger_transcriptome}/genes/genes.gtf.gz ${launchDir}/work
			gunzip ${launchDir}/work/genes.gtf.gz
		fi

		if ! test -f $cellRanger_out/velocyto/*loom; then
  			velocyto run10x -@ `nproc` $cellRanger_out ${launchDir}/work/genes.gtf
		fi
		"""
}

/*
 * produce R markdown output
 * contains initial seurat pipeline
 * can edit .Rmd in output folder and run with -resume to use edited parameters
 */
process seurat_markdown {
	//publishDir 'output', mode: 'copy'
	input:
		val cellranger_out_paths
	output:
		tuple path('**.qs'),path('**markers.txt')
	script:
		"""
  		mkdir -p ${launchDir}/output

		if ! test -f ${launchDir}/output/initial_analysis.Rmd; then
  			cp ${projectDir}/bin/initial_analysis_template.Rmd ${launchDir}/output/initial_analysis.Rmd
		fi

		module load R
		Rscript -e "rmarkdown::render('${launchDir}/output/initial_analysis.Rmd',param=list(args=c('${cellranger_out_paths}')))"
	
		"""
}

/*
 * make shiny app
 */
process make_shiny {
	maxForks 1
	input:
		path seurat_object
	output:
		stdout
	script:
		"""
		if ! test -f ${launchDir}/output/app.R; then
  			cp ${projectDir}/bin/Rshiny_app_template.txt ${launchDir}/output/app.R
		fi

		sh ${projectDir}/bin/create_app_from_template.sh ${launchDir}/output/app.R
		"""
}

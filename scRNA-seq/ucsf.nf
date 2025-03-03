#!/usr/bin/env nextflow

params.expected_cell_number = 10000
params.fastq_source = true
//params.fastq_source = "java -jar ./fdtCommandLine.jar -noupdates -pull -r -c hci-bio-app.hci.utah.edu -d ./ /scratch/fdtswap/fdt_sandbox_gnomex/96ab7c97-8c77-4676-98d8-a43c3bcbf5e7/24250R"
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
		| splitCsv(sep:"\t")
		| set {sample_sheet}

	// get fastq files either from launch folder (default)
	// or from gnomex via fdt
	if(params.fastq_source != true){
		Channel.value(params.fastq_source) 
			| get_fastq
			| flatten 
			| set {fastq_list}
	}
	else {Channel.fromPath("gi*/**fastq.gz") 
		| set {fastq_list}
	}

	// pair up fastq files by sample id
	// replace sample id with name
	fastq_list
		| map{ file ->
            		def sample_id = file.name.replaceAll("_S.*", "")
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
		mkdir -p ${sample_id}_fastqs
		cp $fastqs ${sample_id}_fastqs
		
		#submit to cellranger
		sh ${projectDir}/bin/run_cellranger_ucsf.sh \
			$sample_id \
			${sample_id}_fastqs \
			${params.cellRanger_transcriptome} \
			${params.expected_cell_number}

		#change id to name from here
		mv ${sample_id} ${sample_name}
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

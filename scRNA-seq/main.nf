#!/usr/bin/env nextflow

params.expected_cell_number = 10000
params.fastq_source = true
//params.fastq_source = "empty"
params.sample_table = "*txt"
params.run_velocyto = false
params.miniaturize = false

params.help = false
if (params.help) {
        println """
        scRNA-seq nextflow pipeline

        REQUIRED PARAMETERS:
        --sample_table *txt
                Must be in the launch folder
                Table from gnomex contains "ID" and "Sample Name"
                        gnomex > navigate to experiment > "Experiment Design" > "Download Sample Sheet" > 
                        move to the launch folder
                Can be manually generated with format "ID" "name" "condition"
                        ID is the unique fastq file prefix

        OPTIONAL PARAMETERS:
        --expected_cell_number integer
                Used by CellRanger (default: 10000)

        --run_velocyto true/false
                Run velocyto analysis (default: false)

        --miniaturize true/false
                Create mini fastq files of 250k reads for testing (default: false)

        --fastq_source="source1,source2,..."
                Specify input sources, can combine multiple with commas

                ="java -jar ./fdt....gnomex..."
                        Pulls fastq files from gnomex
                        Get this command from:
                        gnomex > navigate to experiment > "Files" > "Download Files" > 
                        move fastq folder to the right > "FDT Command Line" > copy command

                ="SSD/YYYYMMDD_run_identifier/email_subject_line:password"
                        UCSF core emails a filepath and password after sequencing

                ="SRA"
                        Using sample table text file, downloads fastqs from SRA repository
                        One SRA ID per line (SRR...), skips comments and empty lines

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
            flat_fastq_list | miniaturize | set { fastq_list }
        } else { flat_fastq_list | set { fastq_list } }
	
	// pair up fastq files by sample id
	sample_sheet
		| combine ( fastq_list )
		| filter { sample_id, sample_name, condition, fastq_file ->
			fastq_file.name.startsWith(sample_id + "_")
		}
      		| map { sample_id, sample_name, condition, fastq_file ->
          		def meta = [id: sample_id, name: sample_name, condition: condition]
          		[meta.id, meta, fastq_file]
      		}
      		| groupTuple()  // Group by id
      		| map { id, meta, fastq_files ->
          		[meta[0], fastq_files]
      		}
		| run_cellranger // run pipeline for individual samples

	// Run velocyto splice counting if true
        if (params.run_velocyto == true) {
            run_cellranger.out | run_velocyto | set { alignment_output }
        } else { run_cellranger.out | set { alignment_output } }
	
	// Aggregate and create reports
	alignment_output
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
		Rscript ${projectDir}/../shared_bin/make_sample_table.R $input_sample_table
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
		sh ${projectDir}/../shared_bin/get_files.sh "${input_file_source}"
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
 * run cellRanger
 */
process run_cellranger {
	publishDir 'output', pattern: '**/web_summary.html', mode: 'copy', overwrite: true
	publishDir 'output', pattern: '**/filtered_feature_bc_matrix.h5', mode: 'copy', overwrite: true
	publishDir 'output', pattern: '**/metrics_summary.csv', mode: 'copy', overwrite: true
	input:
		tuple val(meta), path(fastqs, stageAs: 'fastq_folder/*')
	output:
		path "${meta.name}"
	script:
		"""
		module load cellranger
		which cellranger

		cellranger count --id=${meta.name} \
			--fastqs=fastq_folder \
			--transcriptome=${params.cellRanger_transcriptome} \
			--expect-cells=${params.expected_cell_number} \
			--create-bam=${params.run_velocyto} \
			--nosecondary

		#change id to name
		mv ${meta.id} ${meta.name}
		"""
}
		
/*
 * run velocyto
 * produces loom file containing spliced and unspliced information
 */
process run_velocyto {
	publishDir 'output', pattern: '**loom', mode: 'copy', overwrite: true
	input:
		path cellRanger_out
	output:
		path "${cellRanger_out}", emit: file_tree
		path '**loom', emit: loom_file
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

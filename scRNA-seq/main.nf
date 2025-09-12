#!/usr/bin/env nextflow

// Include subworkflows and shared help
include { FASTQ_PREPROCESSING } from '../subworkflows/fastq_preprocessing.nf'
include { WorkflowCompletion } from '../subworkflows/workflow_complete.nf'
include { getSharedHelp } from '../modules/shared_help'

params.expected_cell_number = 10000
params.fastq_source = true
//params.fastq_source = "empty"
params.sample_table = "*txt"
params.run_velocyto = false
params.miniaturize = false
params.outputDir = "output"

params.help = false
if (params.help) {
        println """
        scRNA-seq nextflow pipeline
${getSharedHelp()}

        --expected_cell_number integer
                Used by CellRanger (default: 10000)

        --run_velocyto true/false
                Run velocyto analysis (default: false)

        scRNA-seq specific parameters:
        --adapters.forward, .reverse, .overlap, .nextseqTrim, .minLength, .args
                Adapter trimming parameters (see cutadapt --help)

        """
        exit 0
}

workflow {

	// FASTQ preprocessing (sample table creation, download, decompress, miniaturize, fastq pairing/grouping)
	FASTQ_PREPROCESSING(
		Channel.fromPath(params.sample_table),
		params.fastq_source, 
		params.miniaturize
	)
	
	grouped_fastqs = FASTQ_PREPROCESSING.out.grouped_fastqs
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

workflow.onComplete {
	// Runs on success, cancel or fail
	WorkflowCompletion()
}


/*
 * run cellRanger
 */
process run_cellranger {
	publishDir "${params.outputDir}", pattern: '**web_summary.html', mode: 'copy', overwrite: true
	publishDir "${params.outputDir}", pattern: '**filtered_feature_bc_matrix.h5', mode: 'copy', overwrite: true
	publishDir "${params.outputDir}", pattern: '**metrics_summary.csv', mode: 'copy', overwrite: true
	input:
		tuple val(meta), path(fastqs, stageAs: 'fastq_folder/*')
	output:
		tuple val(meta), path("${meta.name}")
	script:
		"""
		module load cellranger/8.0.1
		
		echo "====== PROCESS_SUMMARY"
		echo "====== RUN_CELLRANGER ======"
		echo "Strategy: Single-cell RNA-seq quantification and analysis"
		echo "cellranger \$(cellranger --version)"
		echo "Expected cells: ${params.expected_cell_number}"
		echo "Transcriptome: ${params.cellRanger_transcriptome}"
		echo "====== RUN_CELLRANGER ======"
		echo "====== PROCESS_SUMMARY"

		echo "== Running cellranger for ${meta.name}"
		cellranger count --id=${meta.name} \\
			--fastqs=fastq_folder \\
			--transcriptome=${params.cellRanger_transcriptome} \\
			--expect-cells=${params.expected_cell_number} \\
			--create-bam=${params.run_velocyto} \\
			--nosecondary > cellranger.log 2>&1
		
		tail -4 cellranger.log
		"""
}
		
/*
 * run velocyto
 * produces loom file containing spliced and unspliced information
 */
process run_velocyto {
	publishDir "${params.outputDir}", pattern: '**loom', mode: 'copy', overwrite: true
	input:
		path cellRanger_out
	output:
		path "${cellRanger_out}", emit: file_tree
		path '**loom', emit: loom_file
	script:
		"""
		module load velocyto
		module load samtools
		
		echo "====== PROCESS_SUMMARY"
		echo "====== RUN_VELOCYTO ======"
		echo "Strategy: Generate spliced/unspliced count matrices for RNA velocity"
		echo "velocyto \$(velocyto --version 2>&1)"
		echo "\$(samtools --version | head -1)"
		echo "====== RUN_VELOCYTO ======"
		echo "====== PROCESS_SUMMARY"

		echo "== Running velocyto for ${cellRanger_out}"
	
		if ! test -f ${launchDir}/work/genes.gtf; then
			echo "== Requires *gtf. Decompressing *gtf.gz into launchDir/work"
			cp ${params.cellRanger_transcriptome}/genes/genes.gtf.gz ${launchDir}/work
			gunzip ${launchDir}/work/genes.gtf.gz
		fi

		if ! test -f $cellRanger_out/velocyto/*.loom; then
  			velocyto run10x -@ ${task.cpus} $cellRanger_out ${launchDir}/work/genes.gtf
			
			if ls $cellRanger_out/velocyto/*.loom 1> /dev/null 2>&1; then
				file_size=\$(du -h $cellRanger_out/velocyto/*.loom | cut -f1)
				loom_file=\$(basename $cellRanger_out/velocyto/*.loom)
				echo "Created \$loom_file: \$file_size"
			else
				echo "== No loom file created"
			fi
		else
			echo "== loom file already exists"
		fi
		"""
}

/*
 * produce R markdown output
 * contains initial seurat pipeline
 * can edit .Rmd in output folder and run with -resume to use edited parameters
 */
process seurat_markdown {
	input:
		val collected_cellranger
	output:
		stdout
	script:
		"""
		echo "====== PROCESS_SUMMARY"
		echo "====== SEURAT_MARKDOWN ======"
		echo "Strategy: Generate Seurat analysis report"
		echo "====== SEURAT_MARKDOWN ======"
		echo "====== PROCESS_SUMMARY"
		
		if ! test -f ${launchDir}/output/initial_analysis.Rmd; then
  			cp ${projectDir}/bin/initial_analysis_template.Rmd ${launchDir}/output/initial_analysis.Rmd
		fi

		module load R
		Rscript -e "rmarkdown::render('${launchDir}/output/initial_analysis.Rmd', \\
			params=list(collected_cellranger='${collected_cellranger}', transcriptome='${params.cellRanger_transcriptome}', velocity='${params.run_velocyto}'))"
	
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
		echo "====== PROCESS_SUMMARY"
		echo "====== MAKE_SHINY ======"
		echo "Strategy: Create interactive Shiny application"
		echo "====== MAKE_SHINY ======"
		echo "====== PROCESS_SUMMARY"
		
		if ! test -f ${launchDir}/output/app.R; then
  			cp ${projectDir}/bin/Rshiny_app_template.txt ${launchDir}/output/app.R
		fi

		sh ${projectDir}/bin/create_app_from_template.sh ${launchDir}/output/app.R
		"""
}

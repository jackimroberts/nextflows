/*
 * FASTQ preprocessing subworkflow
 * Handles sample table creation, file discovery, downloading, decompression, and miniaturization
 */

include { make_sample_table; get_fastq; decompress; miniaturize } from '../modules/shared_processes.nf'

workflow FASTQ_PREPROCESSING {
    take:
    sample_table_file
    params_fastq_source
    params_miniaturize
    
    main:
    // Create sample table
    sample_table_file
        | make_sample_table 
        | splitCsv(sep:"\t") 
        | map { row -> 
            def sample_id = row[0]
            def sample_name = row[1] 
            def condition = row[2]
            def extra_data = row.size() > 3 ? row[3..-1].join('\t') : ""
            [sample_id, sample_name, condition, extra_data]
        }
        | set { sample_sheet }

    // Find existing fastq files
    Channel.fromPath("**fastq*", checkIfExists: false)
        | filter { it.exists() && !it.toString().contains('/work/') } // exclude those in 'work'
        | set { existing_files }

    // Downloaded files if specified and combine with existing files
    // Separate files by type and handle appropriately
    (params_fastq_source != true ? get_fastq(Channel.value(params_fastq_source), sample_table_file) : Channel.empty())
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
        | filter { it.name.matches('.*[-_][RL][12][.-_].*') } // only keep R1/R2/L1/L2 files. I1/I2 are typically the index
        | set { flat_fastq_list }

    // Downsize fastq for faster runs
    if (params_miniaturize == true) {
        flat_fastq_list
            | miniaturize
            | set { fastq_list }
    } else {
        flat_fastq_list
            | set { fastq_list }
    }
    
    emit:
    sample_sheet = sample_sheet
    fastq_files = fastq_list
}
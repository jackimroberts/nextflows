/*
 * FASTQ preprocessing subworkflow
 * Handles sample table creation, file discovery, downloading, decompression, miniaturization, and fastq file pairing/grouping by sample and run
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

    // Pair up fastq files by sample id and run, keeping runs separate
    sample_sheet
        | combine(fastq_list)
        | filter { sample_id, sample_name, condition, extra_data, fastq_file ->
            fastq_file.name.startsWith(sample_id + "_")
        }
        | map { sample_id, sample_name, condition, extra_data, fastq_file ->
            def run_part = fastq_file.name.replaceAll("^${sample_id}_", "").replaceAll("[-_.][LR][12][.-_].*", "")
            def meta = [id: sample_id, name: sample_name, condition: condition, run: run_part, extra: extra_data]
            [[meta.id, meta.run], meta, fastq_file]
        }
        | groupTuple()  // Group by the [id, run] tuple
        | map { id_run_tuple, meta_list, fastq_files ->
            [meta_list[0], fastq_files]
        }
        | set { grouped_fastqs }
    
    emit:
    grouped_fastqs = grouped_fastqs
}
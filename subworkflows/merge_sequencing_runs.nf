/*
 * MERGE_SEQUENCING_RUNS subworkflow
 * 
 * Groups BAM files by sample ID and merges files from multiple sequencing runs
 * while avoiding unnecessary SLURM jobs for single-run samples.
 * 
 * Input: aligned_bams - Channel of [meta, bam_file] from alignment processes
 * Output: merged_bams - Channel of [meta, merged_bam_file] ready for filtering
 */

include { merge_run_bams } from '../modules/shared_processes.nf'

workflow MERGE_SEQUENCING_RUNS {
    take:
    aligned_bams
    
    main:
    // Combine samples by sample id
    aligned_bams
        | map { meta, bam_file ->
            [meta.id, meta, bam_file]
        }
        | groupTuple
        | map { sample_id, meta_list, bam_files ->
            [meta_list[0], bam_files]
        }
    // Separate samples with single or multiple sequencing runs
        | branch {
            multiple: it[1].size() > 1
            single: it[1].size() == 1
        }
        | set { bam_branches }

    // Only merge files from multiple runs
    // Avoids unnecessary SLURM jobs
    bam_branches.multiple 
        | merge_run_bams
    // Combine both streams back together
        | mix( bam_branches.single )
        | set { merged_bams }
    
    emit:
    merged_bams = merged_bams
}
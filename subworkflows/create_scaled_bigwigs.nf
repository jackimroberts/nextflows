/*
 * CREATE_SCALED_BIGWIGS subworkflow
 * 
 * Measures read depth for each BAM file, calculates scaling factors based on minimum depth,
 * and creates normalized BigWig files for visualization.
 * 
 * Input: filtered_bams - Channel of [meta, filtered_bam_file] 
 * Output: bigwigs - Channel of BigWig files for genome browser visualization
 */

include { measure_depth; bam_to_bigwig } from '../modules/shared_processes.nf'

workflow CREATE_SCALED_BIGWIGS {
    take:
    filtered_bams
    
    main:
    // Calculate read depth for each bam
    filtered_bams
        | measure_depth
        | map { meta, bam_file, depth ->
            def depth_value = depth.trim().split('\n').last() as Integer
   
            [meta + [depth: depth_value], bam_file]
        }
        | set { filtered_bams_with_depth }
        
    // Find minimum depth across all samples for scaling
    filtered_bams_with_depth
        | map { meta, bam_file -> meta.depth }
        | min()
    // Combine for bigwig creation
        | combine(filtered_bams_with_depth)
        | bam_to_bigwig
    
    emit:
    bigwigs = bam_to_bigwig.out
}
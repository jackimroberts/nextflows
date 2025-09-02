/*
 * Workflow completion subworkflow
 * Runs resource usage analysis and prints summary
 */

workflow WORKFLOW_COMPLETION {
    take:
    trigger_signal  // Input to trigger completion after main workflow
    
    main:
    // Run completion analysis
    completion_analysis(trigger_signal)
    
    emit:
    completion_report = completion_analysis.out
}

process completion_analysis {
    executor = 'local'
    publishDir "${launchDir}/output", mode: 'copy', overwrite: true
    
    input:
    val trigger_signal
    
    output:
    path "completion_report.txt"
    
    script:
    """
    echo "=========================================================="
    echo "Pipeline execution summary"
    echo "=========================================================="
    echo "Completed at : \$(date)"
    echo "Work directory: ${launchDir}"
    echo "=========================================================="
    
    analyzer_script="${projectDir}/../shared_bin/slurm_usage_analyzer.sh"
    if [ -f "\$analyzer_script" ]; then
        echo "Running resource usage analysis..."
        bash "\$analyzer_script" "${launchDir}" || echo "Resource analysis completed with warnings"
    else
        echo "Resource analyzer script not found at: \$analyzer_script"
        echo "Skipping resource usage analysis."
    fi
    
    echo "=========================================================="
    
    # Save report to file
    {
        echo "Pipeline Completion Report"
        echo "Generated: \$(date)"
        echo "Launch Directory: ${launchDir}"
        echo "Project Directory: ${projectDir}"
        echo ""
        if [ -f "\$analyzer_script" ]; then
            bash "\$analyzer_script" "${launchDir}"
        fi
    } > completion_report.txt
    """
}
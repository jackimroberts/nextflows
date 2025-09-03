// Shared workflow completion handler
// Usage: include { handleWorkflowCompletion } from '../shared_bin/workflow_completion_handler'

def handleWorkflowCompletion() {
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
    
    def analysisOutput = ""
    if (file(analyzerScript).exists()) {
        println "Running resource usage analysis..."
        try {
            def proc = ["bash", analyzerScript, workDir].execute()
            proc.waitFor()
            
            if (proc.exitValue() == 0) {
                analysisOutput = proc.text
                println analysisOutput
            } else {
                analysisOutput = proc.text
                println "Resource analysis completed with warnings:"
                println analysisOutput
                if (proc.err.text) {
                    println "Error output:"
                    println proc.err.text
                    analysisOutput += "\nError output:\n${proc.err.text}"
                }
            }
        } catch (Exception e) {
            analysisOutput = "Failed to run resource analysis: ${e.message}"
            println analysisOutput
        }
    } else {
        analysisOutput = "Resource analyzer script not found at: ${analyzerScript}\nSkipping resource usage analysis."
        println analysisOutput
    }
    
    // Create completion report file
    def outputDir = new File("${workDir}/output")
    if (!outputDir.exists()) {
        outputDir.mkdirs()
    }
    
    def reportFile = new File(outputDir, "completion_report.txt")
    reportFile.text = """Pipeline Completion Report
Generated: ${new Date()}
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
Launch Directory: ${workDir}
Project Directory: ${projectDir}

Resource Usage Analysis:
${analysisOutput}
"""
    println "Completion report saved to: ${reportFile.absolutePath}"
    println "=========================================================="
}
/*
 * Function for handling workflow completion tasks
 * Runs SLURM usage analysis and log collection
 */

def WorkflowCompletion() {
    def launchDir = workflow.launchDir
    def projectDir = workflow.projectDir
    def outputDir = new File("${launchDir}/output")
    if (!outputDir.exists()) outputDir.mkdirs()
    
    def analyzerScript = "${projectDir}/../bin/slurm_usage_analyzer.sh"
    def logScript = "${projectDir}/../bin/collect_task_logs.sh"
    
    def proc1 = ["bash", analyzerScript, launchDir, outputDir.toString()].execute()
    proc1.waitFor()
    
    def proc2 = ["bash", logScript, launchDir, outputDir.toString()].execute()
    proc2.waitFor()
}
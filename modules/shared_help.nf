/*
 * Shared help functions for common parameters
 */

def getSharedHelp() {
    return """
        REQUIRED PARAMETERS:
        --sample_table *txt
                Must be in the launch folder
                Table from gnomex contains "ID" and "Sample Name"
                        gnomex > navigate to experiment > "Experiment Design" > "Download Sample Sheet" > 
                        move to the launch folder
                Can be manually generated with format "ID" "name" "condition"
                        ID is the unique fastq file prefix

        OPTIONAL PARAMETERS:
        --miniaturize true/false
                Create mini fastq files of 2,500 reads for testing (default: false)

        --fastq_source="source1,source2,..."
                Specify input sources, can combine multiple with commas
                
                ="java -jar ./fdt....gnomex..."
                        Pulls fastq files from gnomex
                        Get this command from:
                        gnomex > navigate to experiment > "Files" > "Download Files" > 
                        move fastq folder to the right > "FDT Command Line" > copy command
                
                ="CoreBrowser"
                        Retrieves unarchived files from Utah Core Browser
                        Select files > More options dropdown menu > Secure link >
                        paste into a "core_links" text file
                
                ="SSD/illumina/YYYYMMDD_run_identifier/email_subject_line:password"
                        UCSF core emails a filepath and password after sequencing
                
                ="SRA"
                        Using sample table text file, downloads fastqs from SRA repository
                        Finds any SRR#+, skipping comment and empty lines

                (default: searches for existing fastq files in launch directory)
    """
}
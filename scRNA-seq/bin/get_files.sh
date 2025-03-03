#!/bin/bash
#converts gnomex provided command to utah chpc script
#downloads fastq files

#location of fdt app
FDT="java -jar /uufs/chpc.utah.edu/sys/pkg/fdt/0.9.20/fdt.jar"

#command line argument recieved
file_source=$1

#portion of command to execute
fdt_commands=${file_source##*.jar}
#execute
$FDT$fdt_commands


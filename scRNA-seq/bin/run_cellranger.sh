#!/bin/bash

sample_name=$1
fastqs=$2
transcriptome=$3
cell_number=$4

module load cellranger


/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/10X/cellranger-8.0.1/cellranger count --id=$sample_name \
	--fastqs=$fastqs \
	--transcriptome=$transcriptome \
	--expect-cells=$cell_number \
	--create-bam=false
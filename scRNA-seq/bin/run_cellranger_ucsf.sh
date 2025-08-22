#!/bin/bash

sample_name=$1
fastqs=$2
transcriptome=$3
cell_number=$4

module load cellranger

echo "=== run cellranger count"
which cellranger

cellranger count --id=$sample_name \
	--fastqs=$fastqs \
	--transcriptome=$transcriptome \
	--expect-cells=$cell_number \
	--chemistry=SC3Pv4 \
	--create-bam=true
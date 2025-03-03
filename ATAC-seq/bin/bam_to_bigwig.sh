#!/bin/bash

# convert bam to bigwig
# this requires a scale based on depth of each sample

# rename input variables
sample_num=$1
bam=$2
depth=$3
min_depth=$4

# define application directories
HCI=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl
SAM=$HCI/app/samtools/1.17/samtools
module load deeptools

# create scale based on depth and minimum depth of all samples
scale=$(echo "scale=5; $min_depth/$depth" | bc)

# index bam and generate bigwig
$SAM index $bam
bamCoverage -b $bam --scaleFactor $scale -o $sample_num.bw
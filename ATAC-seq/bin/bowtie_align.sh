#!/bin/bash

FASTQ1=$1
FASTQ2=$2
NAME=$3
INDEX=$4

# load required modules
module load bowtie2
module load samtools

echo "=== align paired-end reads with bowtie2"
which bowtie2
which samtools

bowtie2 --local --very-sensitive --no-mixed --no-discordant -p $(nproc) -x $INDEX \
-1 $FASTQ1 -2 $FASTQ2 | samtools view -bS -o $NAME.raw.bam
#!/bin/bash

HCI=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl
APP=$HCI/app
SAM_APP=$APP/samtools/1.17/samtools

# mouse mm39 index
#INDEX_39=/scratch/general/vast/u0948419/mm39/mm39

# mouse mm10 index
DB=$HCI/data/Mouse/Mm10
INDEX_10=$DB/mm10.standard

FASTQ1=$1
FASTQ2=$2
NAME=$3
module load bowtie2
#which bowtie2
bowtie2 --local --very-sensitive --no-mixed --no-discordant -p `nproc` -x $INDEX_10 \
-1 $FASTQ1 -2 $FASTQ2 | $SAM_APP view -bS -o $NAME.mm10.raw.bam
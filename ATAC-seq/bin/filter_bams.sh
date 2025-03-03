#!/bin/bash

sample_name=$2
filepath=$1
filename=${filepath##*/}
sample=${filename%.raw*}

HCI=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl
CHR=$HCI/data/Mouse/Mm10/mm10.standard.sizes
APP=$HCI/app
SAM=$APP/samtools/1.17/samtools
BED=$HCI/app/bedtools/2.29.0/bedtools
PICARD=$APP/picard/2.26.3/picard.jar
#BLACK=/scratch/general/vast/mm39/mm39.blacklist.bed
BLACK=$HCI/data/Mouse/Mm10/mm10.blacklist.bed

$SAM sort -m 4G $filepath -o $sample.sorted.bam

#mark and remove duplicates
java -jar $PICARD MarkDuplicates \
  -I $sample.sorted.bam -O $sample.dedup.bam -M $sample.dup_metrics.txt \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --TAGGING_POLICY All \
  --REMOVE_DUPLICATES \
  --VERBOSITY WARNING

# filter out any duplicate, unaligned, secondary alignments
#I've chosen to leave reads without mates because it's atac-seq
#Normally it would be -F 3340
$SAM view -F 3332 $sample.dedup.bam -o $sample.filt.bam

#remove mitochondrial reads
$SAM index $sample.filt.bam
cat $CHR | cut -f 1 | grep -v M | xargs $SAM view -b $sample.filt.bam > $sample.tmp

#remove blacklisted regions
$BED intersect -v -abam $sample.tmp -b $BLACK > $sample.filtered.bam

#calculate 
echo $($SAM view -f 64 -c $sample.filtered.bam) | tr -d '\n'
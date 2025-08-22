#!/bin/bash

sample_name=$2
filepath=$1
filename=${filepath##*/}
sample=${filename%.raw*}

CHR=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Mouse/Mm10/mm10.standard.sizes
#BLACK=/scratch/general/vast/mm39/mm39.blacklist.bed
BLACK=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/Mouse/Mm10/mm10.blacklist.bed
PICARD=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/picard/2.26.3/picard.jar

module load samtools
module load bedtools

echo "=== filter bams for duplicate, unaligned, secondary, mitochondrial or blacklisted"
which samtools
which bedtools

samtools sort -m 4G $filepath -o $sample.sorted.bam

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
samtools view -F 3332 $sample.dedup.bam -o $sample.filt.bam

#remove mitochondrial reads
samtools index $sample.filt.bam
cat $CHR | cut -f 1 | grep -v M | xargs samtools view -b $sample.filt.bam > $sample.tmp

#remove blacklisted regions
bedtools intersect -v -abam $sample.tmp -b $BLACK > $sample.filtered.bam

#calculate 
echo $(samtools view -f 64 -c $sample.filtered.bam) | tr -d '\n'
#!/bin/bash

sample=$2
filepath=$1
blacklist=$3

module load samtools
module load bedtools

echo "=== filter bams for duplicate, unaligned, secondary, mitochondrial or blacklisted"
which samtools
which bedtools

# Sort BAM file by name for fixmate
samtools sort -n -@ $(nproc) -m 4G $filepath -o $filepath.tmp

# Run fixmate to fill in mate coordinates and insert size fields
samtools fixmate -m $filepath.tmp $filepath.fixmate.bam

# Sort by coordinate for markdup
samtools sort -@ $(nproc) -m 4G $filepath.fixmate.bam -o $filepath.sorted.bam

# Mark duplicates including optical duplicates using samtools
samtools markdup -d 2500 -@ $(nproc) $filepath.sorted.bam $filepath.dedup.bam

# Index the deduplicated BAM file for idxstats
samtools index $filepath.dedup.bam

# Get non-mitochondrial chromosome list
CHROMS=$(samtools idxstats $filepath.dedup.bam | cut -f 1 | grep -v -E '(MT|chrM|chrMT)' | tr '\n' ' ')

# Filter out duplicates, unaligned, secondary alignments, and chrM in one step
samtools view -@ $(nproc) -b -F 3332 $filepath.dedup.bam $CHROMS > $filepath.filt.bam

# Index intermediate file
samtools index $filepath.filt.bam

# Remove blacklist regions if provided
if [ -f "$blacklist" ]; then
	bedtools intersect -v -abam $filepath.filt.bam -b $blacklist > $filepath.filtered.bam
else
    cp $filepath.filt.bam $filepath.filtered.bam
fi

# Index final filtered BAM
samtools index $filepath.filtered.bam

# Clean up intermediate files
rm -f $filepath.tmp $filepath.fixmate.bam $filepath.sorted.bam

# Calculate size for eventual scaling
echo $(samtools view -f 64 -c $filepath.filtered.bam) | tr -d '\n'

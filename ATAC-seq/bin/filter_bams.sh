#!/bin/bash

## ---------------------------
##
## Script name: filter_bams
##
## Purpose of script: 
##	Filter BAM files for duplicates, unaligned reads, secondary alignments,
##	mitochondrial reads, and blacklisted regions. Returns final read count.
##
## ---------------------------
##
## Usage:
##	./filter_bams.sh <filepath> <sample_id> <blacklist>
##	
##	filepath: Path to input BAM file
##	sample_id: Sample identifier
##	blacklist: Path to blacklist BED file (optional)
##
## ---------------------------
##
## Notes:
##	Uses samtools for duplicate marking and filtering
##	Uses bedtools for blacklist region filtering
##	Optimized for ATAC-seq data processing
##
## ---------------------------

## Process command line arguments
filepath=$1
sample=$2
blacklist=$3

## Load required modules
module load samtools
module load bedtools

echo "====== PROCESS_SUMMARY"
echo "====== FILTER_BAMS ======"
echo "Strategy: Filter BAM files"
echo "mark optical duplicates with -d 2500"
echo "remove duplicate, secondary and unaligned with -F 3332"
echo "remove mitochondrial reads"
echo "samtools $(samtools --version | head -1)"
echo "bedtools $(bedtools --version)"
echo "blacklist: $blacklist"
echo "====== FILTER_BAMS ======"
echo "====== PROCESS_SUMMARY"

echo "filtering: $filepath"

## Count initial reads
echo "=== Initial Read Counts"
initial_reads=$(samtools view -c $filepath)
echo "Initial reads: $initial_reads"

## Sort BAM file by name for fixmate
samtools sort -n -@ $(nproc) $filepath -o $sample.tmp

## Run fixmate to fill in mate coordinates and insert size fields
samtools fixmate -m $sample.tmp $sample.fixmate.bam

## Sort by coordinate for markdup
samtools sort -@ $(nproc) $sample.fixmate.bam -o $sample.sorted.bam

## Mark duplicates including optical duplicates using samtools
samtools markdup -d 2500 -@ $(nproc) $sample.sorted.bam $sample.dedup.bam

## Count duplicates
echo "=== Duplicate Statistics"
total_after_dedup=$(samtools view -c $sample.dedup.bam)
duplicates=$(($initial_reads - $total_after_dedup))
echo "Duplicates marked: $duplicates"
echo "Reads after deduplication: $total_after_dedup"

## Index the deduplicated BAM file for idxstats
samtools index $sample.dedup.bam

## Get non-mitochondrial chromosome list
CHROMS=$(samtools idxstats $sample.dedup.bam | cut -f 1 | grep -v -E '(MT|chrM|chrMT)' | tr '\n' ' ')

## Count mitochondrial reads before filtering
echo "=== Mitochondrial Read Statistics ==="
mito_reads=$(samtools idxstats $sample.dedup.bam | grep -E '(MT|chrM|chrMT)' | awk '{sum+=$3} END {print sum+0}')
echo "Mitochondrial reads: $mito_reads"

## Filter out duplicates, unaligned, secondary alignments, and chrM
samtools view -@ $(nproc) -b -F 3332 $sample.dedup.bam $CHROMS > $sample.filt.bam

## Count reads after filtering
echo "=== Post-filtering Statistics"
reads_after_filter=$(samtools view -c $sample.filt.bam)
echo "Reads after removing duplicates, unaligned, secondary, and mitochondrial: $reads_after_filter"

## Index intermediate file
samtools index $sample.filt.bam

## Remove blacklist regions if provided
if [ -f "$blacklist" ]; then
	echo "=== Blacklist Filtering"
	bedtools intersect -v -abam $sample.filt.bam -b $blacklist > $sample.filtered.bam
	blacklist_removed=$(($reads_after_filter - $(samtools view -c $sample.filtered.bam)))
	echo "Reads removed by blacklist filtering: $blacklist_removed"
else
	echo "=== No Blacklist Filtering"
	echo "No blacklist file provided, skipping blacklist filtering"
    cp $sample.filt.bam $sample.filtered.bam
fi

## Index final filtered BAM
samtools index $sample.filtered.bam

## Final read count summary
echo "=== Final Summary"
final_reads=$(samtools view -c $sample.filtered.bam)
echo "Final filtered reads: $final_reads"
echo "Total reads removed: $(($initial_reads - $final_reads))"
echo "Retention rate: $(echo "scale=2; $final_reads * 100 / $initial_reads" | bc)%"

## Clean up intermediate files
rm -f $sample.tmp $sample.fixmate.bam $sample.sorted.bam
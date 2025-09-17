g#!/bin/bash

## ---------------------------
##
## Script name: filter_bams
##
## Purpose of script: 
##	Filter BAM files for duplicates, unaligned reads, secondary alignments,
##	unmated reads (optional), mitochondrial reads, and blacklisted regions.
##	Filtering controlled by samtools flag parameter.
##
## ---------------------------
##
## Usage:
##	./filter_bams.sh <filepath> <sample_id> <blacklist> <flag> <cpus>
##	
##	filepath: Path to input BAM file
##	sample_id: Sample identifier
##	blacklist: Path to blacklist BED file (optional)
##	flag: Samtools flag value for filtering (e.g., 3340, 3332)
##	cpus: Number of CPU threads
##
## ---------------------------
##
## Notes:
##	Uses samtools for duplicate marking and filtering
##	Uses bedtools for blacklist region filtering
##
## ---------------------------

## Process command line arguments
filepath=$1
sample=$2
blacklist=$3
flag=$4
cpus=$5
arg_markdup=$6
## see flags: https://broadinstitute.github.io/picard/explain-flags.html

## Load required modules
module load samtools
module load bedtools

echo "====== PROCESS_SUMMARY"
echo "====== FILTER_BAMS ======"
echo "Strategy: Filter BAM files"
echo "mark optical duplicates with -d 2500"
if [ -n "$arg_markdup" ]; then
    echo "markdup args: $arg_markdup"
fi
echo "filter with -F $flag"
if [ $flag == 3340 ]; then
    echo "remove duplicate, secondary, unaligned and unmated"
elif [ $flag == 3332 ]; then
    echo "remove duplicate, secondary and unaligned"
fi
echo "remove mitochondrial reads"
if [ -f "$blacklist" ]; then
    echo "blacklist: $blacklist"
fi
echo "samtools $(samtools --version | head -1)"
echo "bedtools $(bedtools --version)"
echo "====== FILTER_BAMS ======"
echo "====== PROCESS_SUMMARY"

echo "filtering: $filepath"

## Sort BAM file by name
## Run fixmate to fill in mate coordinates and insert size fields
## Sort by coordinate
## Mark duplicates including optical duplicates
samtools sort -n -@ $cpus $filepath | \
samtools fixmate -@ $cpus -m - - | \
samtools sort -@ $cpus - | \
samtools markdup $arg_markdup -s -d 2500 -@ $cpus - $sample.dup.bam 2> $sample.duplicate.stats

samtools index $sample.dup.bam
samtools idxstats $sample.dup.bam | sort -V > $sample.idxstats

## Count reads efficiently in single pass
read_counts=$(samtools view $sample.dup.bam | awk '
BEGIN {total=0; unaligned=0; secondary=0; unmated=0}
{
    total++
    flag = $2
    if (and(flag, 4)) unaligned++
    if (and(flag, 256)) secondary++
    if (and(flag, 8)) unmated++
}
END {
    print total " " unaligned " " secondary " " unmated
}')

# Extract individual counts and display
initial_reads=$(echo $read_counts | cut -d' ' -f1)
unaligned_reads=$(echo $read_counts | cut -d' ' -f2)
secondary_reads=$(echo $read_counts | cut -d' ' -f3)
unmated_reads=$(echo $read_counts | cut -d' ' -f4)

echo "$initial_reads: Initial reads"
echo "$unaligned_reads : unaligned"
echo "$secondary_reads : secondary"
echo "$unmated_reads : unmated"

## Parse markdup statistics from stderr output
total_duplicates=$(grep "DUPLICATE TOTAL:" $sample.duplicate.stats | awk '{print $3}')
optical_paired=$(grep "DUPLICATE PAIR OPTICAL:" $sample.duplicate.stats | awk '{print $4}')
optical_single=$(grep "DUPLICATE SINGLE OPTICAL:" $sample.duplicate.stats | awk '{print $4}')
optical_duplicates=$((optical_paired + optical_single))
other_duplicates=$((total_duplicates - optical_duplicates))

echo "$optical_duplicates : optical duplicates (${optical_paired} paired + ${optical_single} single)"
echo "$other_duplicates : other duplicates"

## Get non-mitochondrial chromosome list from stored idxstats
CHROMS=$(cat $sample.idxstats | cut -f 1 | grep -v -E '(MT|chrM|chrMT)' | tr '\n' ' ')

## Count mitochondrial reads from stored idxstats
mito_reads=$(cat $sample.idxstats | grep -E '(MT|chrM|chrMT)' | awk '{sum+=$3} END {print sum+0}')
echo "$mito_reads : mitochondrial"

## Filter out duplicates, unaligned, secondary alignments, and chrM
samtools view -@ $cpus -b -F $flag $sample.dup.bam $CHROMS > $sample.filt.bam

## Count reads after filtering
reads_after_filter=$(samtools view -c $sample.filt.bam)

## Index intermediate file
samtools index $sample.filt.bam

## Remove blacklist regions if provided
if [ -f "$blacklist" ]; then
	# Check if chromosome names match between BAM and blacklist
	bam_chr=$(head -1 $sample.idxstats | cut -f1)
	blacklist_chr=$(head -1 $blacklist | cut -f1)

	if [[ "$bam_chr" == chr* && "$blacklist_chr" != chr* ]]; then
		# Add chr prefix to blacklist
		awk '{print "chr"$0}' $blacklist > blacklist.fixed.bed
		blacklist_file=blacklist.fixed.bed
	elif [[ "$bam_chr" != chr* && "$blacklist_chr" == chr* ]]; then
		# Remove chr prefix from blacklist
		sed 's/^chr//' $blacklist > blacklist.fixed.bed
		blacklist_file=blacklist.fixed.bed
	else
		blacklist_file=$blacklist
	fi

	bedtools intersect -v -abam $sample.filt.bam -b $blacklist_file > $sample.filtered.bam
	blacklist_removed=$(($reads_after_filter - $(samtools view -c $sample.filtered.bam)))
	echo "$blacklist_removed : blacklisted"
else
	echo "No blacklist file provided, skipping blacklist filtering"
    cp $sample.filt.bam $sample.filtered.bam
fi

## Index final filtered BAM
samtools index $sample.filtered.bam

## Final read count summary
final_reads=$(samtools view -c $sample.filtered.bam)
echo "Reads retained: $final_reads"
echo "Total reads removed: $(($initial_reads - $final_reads))"
echo "Retention rate: $(echo "scale=2; $final_reads * 100 / $initial_reads" | bc)%"
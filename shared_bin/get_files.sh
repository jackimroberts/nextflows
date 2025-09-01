#!/bin/bash

## ---------------------------
##
## Script name: get_files
##
## Purpose of script: 
##	Download fastq files from multiple sources (GNomEx, Core Browser, UCSF, SRA)
##	Supports parallel downloads and generates checksums
##
## ---------------------------
##
## Usage:
##	./get_files.sh <input_sources> <sample_table> <launch_dir>
##	
##	input_sources: Comma-separated list of sources (gnomex, CoreBrowser, UCSF:password, SRA)
##	sample_table: Path to sample table file (required for SRA downloads)
##	launch_dir: Launch directory path (for locating core_links file)
##
## ---------------------------
##
## Notes:
##	Handles comma-separated list of sources
##	Uses appropriate tools for each source: FDT (GNomEx), aria2 (Core Browser), lftp (UCSF), sra-toolkit (SRA)
##	Optimized for parallel processing using available CPU cores
##
## ---------------------------


## Handle command line arguments received
input_sources="$1"
sample_table="$2"
launch_dir="$3"

## Load required modules
module load parallel

## Split comma-separated sources and process each
IFS=',' read -ra SOURCES <<< "$input_sources"
for source in "${SOURCES[@]}"; do
	# Trim whitespace
      	source=$(echo "$source" | xargs)
      	echo "Processing source: $source"


	if [[ "$source" == *"gnomex"* ]]; then
		echo "=== download fastq files from GNomEx, using FDT command line"

		# true location of fdt app with parallel optimization
		FDT="java -jar /uufs/chpc.utah.edu/sys/pkg/fdt/0.9.20/fdt.jar"

		# portion of command to execute
		fdt_commands="${source##*.jar}"
		# execute with parallel streams
		$FDT -P $(nproc) $fdt_commands

	elif [[ "$source" == "CoreBrowser" ]] && [[ -f "${launch_dir}/core_links" ]]; then
		echo "=== download fastq files from Utah core browser via aria"

		module load aria2

		aria2c -i "${launch_dir}/core_links" -j $(nproc) -x 16

	elif [[ "$source" == *":"* ]]; then
		echo "=== download fastq files from UCSF core"
	
		# parse directory and password from format: directory:password
		directory="${source%%:*}"
		password="${source##*:}"

		# use lftp for parallel transfer
		lftp -u hiseq_user,"$password" sftp://fastq.ucsf.edu <<-EOF
			set cmd:parallel $(nproc)
			set net:max-retries 3
			cd $directory
			mget *
			quit
		EOF

	elif [[ "$source" == "SRA" ]] && [[ -f "$sample_table" ]]; then
		echo "=== download fastq files from SRA"
	
		# load SRA toolkit
		module load sra-toolkit
	
		# Extract SRA IDs (skip comments/blank lines, look for SRR pattern)
      		grep -v '^#' "$sample_table" | grep -v '^[[:space:]]*$' | \
      		grep -oE 'SRR[0-9]+' | \
		parallel -j$(nproc) \
			'echo "Downloading {}..."; prefetch {} && fasterq-dump {} --split-3 --gzip --skip-technical && \
				echo "Completed {}" || echo "FAILED: {}"'
	
	else
		echo "=== couldn't identify source"
	fi
done

## Generate checksums for all downloaded files
echo "=== Generating checksums"
md5sum **fastq* > md5_downloads.txt 2>/dev/null || echo "No fastq files to checksum"
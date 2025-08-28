#!/bin/bash
# converts gnomex provided command to utah chpc script
# downloads fastq files

# command line argument received
file_source="$1"

module load parallel

if [[ "$file_source" == *"gnomex"* ]]; then
	echo "=== download fastq files from GNomEx, using FDT command line"

	# true location of fdt app with parallel optimization
	FDT="java -jar /uufs/chpc.utah.edu/sys/pkg/fdt/0.9.20/fdt.jar"

	# portion of command to execute
	fdt_commands="${file_source##*.jar}"
	# execute with parallel streams
	$FDT -P $(nproc) $fdt_commands

elif [[ "$file_source" == *":"* ]]; then
	echo "=== download fastq files from UCSF core"
	
	# parse directory and password from format: directory:password
	directory="${file_source%%:*}"
	password="${file_source##*:}"

	# use lftp for parallel transfer
	lftp -u hiseq_user,"$password" sftp://fastq.ucsf.edu <<-EOF
		set cmd:parallel $(nproc)
		set net:max-retries 3
		cd $directory
		mget *
		quit
	EOF

elif [[ "$file_source" == "SRA" ]] && [[ -f "SRA_ids" ]]; then
	echo "=== download fastq files from SRA"
	
	# load SRA toolkit
	module load sra-toolkit
	
	# Filter out comments and empty lines, then download in parallel
	grep -v '^#' SRA_ids | grep -v '^$' | \
	parallel -j$(nproc) \
		'echo "Downloading {}..."; prefetch {} && fasterq-dump {} --split-3 --gzip --skip-technical && echo "Completed {}" || echo "FAILED: {}"'
	
else
	echo "=== couldn't identify source"
	echo "Expected formats:"
	echo "  GNomEx: 'java -jar ./fdt...gnomex...'"
	echo "  UCSF:   'directory_path:password'"
	echo "  SRA:    'SRA' (requires SRA_ids file in current directory)"
fi

# Generate checksums for all downloaded files
echo "=== Generating checksums"
md5sum **/*fastq* > md5.txt 2>/dev/null || echo "No fastq files to checksum"
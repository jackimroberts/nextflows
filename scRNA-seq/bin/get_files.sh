#!/bin/bash
# converts gnomex provided command to utah chpc script
# downloads fastq files

# command line argument received
file_source="$1"

if [[ "$file_source" == *"gnomex"* ]]; then
	echo "=== download fastq files from GNomEx, using FDT command line"

	# true location of fdt app
	FDT="java -jar /uufs/chpc.utah.edu/sys/pkg/fdt/0.9.20/fdt.jar"

	# portion of command to execute
	fdt_commands="${file_source##*.jar}"
	# execute
	$FDT $fdt_commands

elif [[ "$file_source" == *":"* ]]; then
	echo "=== download fastq files from UCSF core"
	
	# parse directory and password from format: directory:password
	directory="${file_source%%:*}"
	password="${file_source##*:}"

	# use lftp for transfer
	lftp -u hiseq_user,"$password" sftp://fastq.ucsf.edu <<-EOF
		cd $directory
		mget *
		quit
	EOF

	# Generate checksums
	md5sum * > md5.txt 2>/dev/null || echo "No files to checksum"
else
	echo "=== couldn't identify source"
	echo "Expected formats:"
	echo "  GNomEx: 'java -jar ./fdt...gnomex...'"
	echo "  UCSF:   'directory_path:password'"
fi
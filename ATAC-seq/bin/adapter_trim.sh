#!/bin/bash

HCI=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl

# adapters for standard compatible Illumina TruSeq
ADAPTF=CTGTCTCTTATACACATCT
ADAPTR=CTGTCTCTTATACACATCT

# application paths
APP=$HCI/app
CUTADAPT=$APP/modulesoftware/cutadapt

#cutadapt
echo "=== adapter trimming"
which $CUTADAPT
$CUTADAPT -O 1 --nextseq-trim=20 -m 1 -a $ADAPTF \
	-A $ADAPTR \
	-o $1.1.fq -p $1.2.fq $2 $3

#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: make_sample_table
##
## Purpose of script: 
##	generate consistent sample detail table for pipeline
##
## ---------------------------
##
## Usage:
##	./make_sample_table.R <input_filepath>
##	
##	input_filename: Path to input sample table file
##
## ---------------------------
##
## Notes:
##	Tested for sample table from Gnomex
##	Also works for manually curated table with ID, name, condition
##
## ---------------------------

## load packages:

require(tidyverse)

## ---------------------------

## manipulate table to have c(ID,Sample_Name,condition,...)
## ID should match the fastq and be unique
## Sample_Name and condition will be generated if not present
##

args = commandArgs(TRUE)

input_filepath <- args[1]

# read input_table
# exclude commented lines, assume header
input_table<-read.delim(input_filepath,
	comment.char="#",
	sep="\t")

# determine if the colnames is likely a row
has_header = FALSE 
for (test_col in 1:ncol(input_table)){
	# if all values in a column are the length, and the column name is different
	# then it's probably header
	if (length (unique (str_length (input_table[,test_col])))==1 
		& !identical(str_length(input_table[1,test_col]),str_length(colnames(input_table)[test_col]))) {
		has_header = TRUE
	}
}

# move colnames into the table if there's no header
if (has_header==FALSE){ input_table <- rbind(colnames(input_table),input_table) }

# assume there is always an id for fastq files first
# if there is only id, use that for name as well
if (ncol(input_table) == 1){ input_table[,2] <- input_table[,1] }

# make sure name doesn't contain spaces
input_table[,2] <- str_replace_all(input_table[,2]," ","_")

# if there is no condition, use name
# take only information after a delimiter.
if (ncol(input_table) == 2){ input_table[,3] <- sub(".*[_-]","",input_table[,2]) }

## ---------------------------

## write table

write_tsv(input_table,"sample_table.tsv",col_names=FALSE)

input_table
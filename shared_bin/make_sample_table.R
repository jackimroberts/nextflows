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
##	./make_sample_table.R <input_filename>
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

## manipulate table to have c(ID,Sample_Name,condition)
## ID should match the fastq and be unique
## Sample_Name is arbitrary and can be the same as ID
## condition will be compared

args = commandArgs(TRUE)

input_filename <- args[1]

input_table <- read_tsv(input_filename) 

if(all(c("ID","Sample Name") %in% colnames(input_table))){
	input_table<- input_table %>% 
		select(ID,`Sample Name`) %>%
		rename(name=`Sample Name`) %>%
		mutate(name=str_replace_all(name," ","_"), # spaces are now underscores
			name=str_replace_all(name,"#",""), # removed hashtags
			condition = sub(".*_","",name))	# took the last part of the name as condition...
}

## ---------------------------

## write table

write_tsv(input_table,"sample_table.tsv",col_names=FALSE)

print(input_table)
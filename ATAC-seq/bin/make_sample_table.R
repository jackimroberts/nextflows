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
## Notes:
##	Currently only tested for sample table from Gnomex
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

input_table <- read_tsv(input_filename) %>% 
	select(ID,`Sample Name`) %>%
	rename(name=`Sample Name`) %>%
	mutate(name=str_replace_all(name," ","_"), # spaces are now underscores
		name=str_replace_all(name,"#",""), # removed hashtags
		condition = sub(".* ","",name))	# took the last part of the name as condition...

## ---------------------------

## write table

write_tsv(input_table,"sample_table.tsv",col_names=FALSE)
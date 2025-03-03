#!/usr/bin/bash

printf "\n#don't forget to\nmodule load nextflow\n\n"
printf "\n#I recommend you run everything in an interactive slurm session\nor you can set up the config file to run on slurm\n\n"
printf '\n#run nextflow pipeline with:\nnextflow *nf\n#or when directing to a folder:\nnextflow run folder\n\n'

printf '\n#get details and parameters with: \nnextflow *nf --help\n\n'

printf "\n#to avoid repeating steps: \nnextflow *nf -resume\n\n"

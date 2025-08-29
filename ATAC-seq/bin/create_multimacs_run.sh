#!/bin/bash

template_file=$1

chip_details=$(cat - | tr "]" " " | tr "[" " " | tr "\n" " " | sed -e 's/, --/--/')

cat $template_file \
	| sed -e "s~#details_here~$chip_details \\\~" > multimacs_command.sh
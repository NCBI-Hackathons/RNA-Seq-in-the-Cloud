#!/bin/bash

for file in `ls -1Sr tight_novel_exon_novel*.gtf.tsv.tsv`; do
	cnt=$(grep -Poc 'STRG\.\d+\.\d+' $file) 
	echo $file $cnt
done

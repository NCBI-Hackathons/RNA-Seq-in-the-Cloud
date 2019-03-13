#!/bin/bash
for file in `ls -1Sr exon_novel*.gtf.tsv`;do
	echo $file
	python TSS_support.py $file CAGE_clusters_C.gff3
	echo 'done this one'
done

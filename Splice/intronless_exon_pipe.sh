#!/bin/bash

for file in *.gtf; do
	python intronless_exon.py $file
done

#for file in `ls -lSr exon_novel_*.gtf.tsv`; do
#	python TSS_support.py $file CAGE_clusters_C.gff3
#done

#!/bin/bash

for f in *.merged.bam ; do
    study=$(echo $f | cut -f1 -d'_')
    tissue=$(echo $f | cut -f2 -d'_' | sed 's/\.merged\.bam//g')
    track=$(echo $f | sed 's/\.merged\.bam//g')
    label="${study}, ${tissue}"
    bigDataUrl="https://storage.googleapis.com/ncbi_sra_rnaseq/${f}"
    echo "track $track"
    echo "bigDataUrl $bigDataUrl"
    echo "shortLabel $label"
    echo "longLabel $label"
    echo "type bam"
done 

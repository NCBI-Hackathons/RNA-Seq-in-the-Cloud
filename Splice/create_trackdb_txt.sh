#!/bin/bash

base=$1

for f in ${base} ; do
    study=$(echo $f | sed -r 's/([^_]*).*/\1/g')
    tissue=$(echo $f | sed -r 's/([^_]*)_(.*)\.merged\.bam/\2/g')
    track=$(echo $f | sed 's/\.merged\.bam//g')
    label="${study}, ${tissue}"
    bigDataUrl="https://storage.googleapis.com/ncbi_sra_rnaseq/${f}"
    echo "track $track"
    echo "bigDataUrl $bigDataUrl"
    echo "shortLabel $label"
    echo "longLabel $label"
    echo "type bam"
done

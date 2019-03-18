#!/bin/bash -eu
# helper to merge transcripts blast outputs into a common cvs file

find /data/contigs/blast_transcripts/ -name '*.transcripts' | while read t; do acc=$(echo $t | tr '[/.]' '\t' | cut -f 5); grep -v '^#' $t | sed "s/^/$acc\t/" ; done | tr '\t' ',' > transcripts.csv

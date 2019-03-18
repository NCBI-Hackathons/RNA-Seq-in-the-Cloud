#!/bin/sh -eu
# helper to merge lnc blast outputs into a common cvs file

find /data/contigs/blast_lnc/ -name '*.lnc' | while read t; do acc=$(echo $t | tr '[/.]' '\t' | cut -f 5); grep -v '^#' $t | sed "s/^/$acc\t/" ; done | tr '\t' ',' > lnc.cvs

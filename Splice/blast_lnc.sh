#!/bin/bash -eu
# blast fasta file against long non-coding sequences from lncipedia 5.2 and capture the output table into *.lnc file

fasta=$1
blastn -query $fasta -db /opt/refs/lncipedia_5_2.fasta -outfmt '7 qseqid sseqid sacc slen length sstart send pident evalue' -num_alignments 5  -perc_identity 95 > $(basename $fasta).lnc

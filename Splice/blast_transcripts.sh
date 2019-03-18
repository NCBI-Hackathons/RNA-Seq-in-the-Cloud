#!/bin/bash -eu
# blast fasta file againasd human trascripts database and capture output table into *.transcripts file

fasta=$1
blastn -query $fasta -db /opt/refs/human_transcripts.fa -outfmt '7 qseqid sseqid sacc slen length sstart send pident evalue' -num_alignments 5  -perc_identity 95 > $(basename $fasta).transcripts

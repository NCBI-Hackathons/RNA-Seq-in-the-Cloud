#!/bin/bash -eu
# calculate sequence md5 for a fasta file and output to *.md5 file

fasta=$1
acc=$(echo $fasta | sed -r 's/.*([ESD]RR[0-9]+).*/\1/')
bioawk -cfastx '{print $name, length($seq), $seq}' $fasta | while read -r name len seq; do
  md5=$(echo -n "$seq" | md5sum - | cut -f 1 -d ' ')
  echo -e "$acc\t$name\t$len\t$md5"
  done > $acc.md5

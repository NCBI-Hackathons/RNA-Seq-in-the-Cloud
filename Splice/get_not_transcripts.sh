#!/bin/bash -eu
# given a contig fasta it matches contigs that are mentioned in the corresponding transcripts file
# and prints only those contigs that are not mentioned in the blast output

fasta=$1
transcripts=$(echo $fasta | sed -r 's/.*\/([ESD]RR.*)$/\1.transcripts/')
diff --new-line-format="" --unchanged-line-format=""  <(cut -f 1 $fasta.fai | sort) <(grep -v '^#' $transcripts | cut -f1 | sort)
exit 0

find /data/contigs/ -maxdepth 1 -name '*.contigs.fa' | while read fasta; do 
  export trasnscripts=$(echo $fasta | sed -r 's/.*\/([ESD]RR.*)$/\1.transcripts/')
  cut -f 1 ${fasta}.fai | while read contig; do
    if ! grep -q $contig $trasnscripts; then echo QQQ; fi
  done
done

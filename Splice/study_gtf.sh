#!/bin/bash -eux

gtfin=$1
threads=${2:-4}
while read line; do
  study=$(echo "$line" | cut -f1)
  tissue=$(echo "$line" | cut -f2 | tr ' ' '_' )
  runs=$(echo "$line" | cut -f3)
  mkdir $study>/dev/null ||:
  cd $study
  gsutil -m cp -P -n $(for run in $(echo "$runs" | tr , ' '); do
    echo -ne " gs://ncbi_sra_rnaseq/$run.bam"
  done) . 
  tmpbam=$(mktemp merged.XXX.bam)
  samtools merge -@4 -f $tmpbam *.bam
  tmpbam='merged.84a.bam'
  samtools index $tmpbam
  chrbam=$(mktemp chr.XXX.bam)
  samtools view -h $tmpbam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | samtools view -@$threads -Sb - > $chrbam
  gtfout="${study}_${tissue}.gtf"
  stringtie $chrbam -p $threads -G $gtfin -o $gtfout
done


#!/bin/bash -eux

function build_gtf(){
  local study=$1
  shift
  local tissue=$1
  shift
  local threads=$1
  shift
  local runs=$@

  if ! gsutil -m cp -P -n $(for run in $(echo "$runs" | tr , ' '); do
    echo -ne " gs://ncbi_sra_rnaseq/$run.bam"
  done) . ; then
    if ! ls *.bam; then
      echo "Failed to pull any BAM for $study $tissue" >&2
      return 1
    fi
  fi
  for bam in *[0-9].bam; do
    chrbam=${bam/%.bam/.chr.bam}
    if [ -f $chrbam ]; then continue; fi
    samtools index $bam
    samtools view -h $bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | samtools view -@$threads -Sb - > $chrbam
#    rm $bam
  done
  mergedbam="${study}_${tissue}.bam"
  samtools merge -@4 -f - *.chr.bam | samtools sort -@8 - > $mergedbam
  samtools index $mergedbam

  gtfout="${study}_${tissue}.gtf"
  stringtie $mergedbam -p $threads -G $gtfin -o $gtfout
  echo $mergedbam $gtfout
}

gtfin=$1
threads=${2:-4}
while read line; do
  study=$(echo "$line" | cut -f1)
  tissue=$(echo "$line" | cut -f2 | tr ' ' '_' )
  runs=$(echo "$line" | cut -f3)
  mkdir $study>/dev/null ||:
  pushd $study
  
  build_gtf $study $tissue $threads $runs ||:

  popd
  exit 1
done


#Script by: SpliceItUp group
#UNC - Chapel Hill 
#Hackathon: RNA-seq-in-the-cloud

#Running this script locally:   python3 analyze_intron_structure.py

import csv
import numpy as np
import pandas as pd

import os
import sys
import gzip
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL

#path_to_files = "/home/vkkodali/example-gtf-files/"
#path_to_files = "/home/ijimenez/RNA-Seq-in-the-Cloud/Splice/"
path_to_files = "/Users/ivanjimenezruiz/Desktop/CLASES/UNC/Hackathons/RNA_Seq_in_the_cloud/"
names=[]
ages=[]
rejectTx = []
intronless_Tx = []

def write_gtf_from_array(list_txs,file_name):
    with open(file_name, 'w') as file:
        for item in list_txs:
            item = '\t'.join(item)
            file.write("%s\n" % item)

with open(path_to_files + 'SRR1045766.bam.GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.harmonized_utrs.gtf','r') as bamfile:

    #with open(path_to_files + 'a_subset.gtf','r') as bamfile:
    reader = csv.reader(bamfile, delimiter = '\t')
    
    # check the header
    finding_intron = 0
    exon_counts = 0
    previous_tx = ""
    rejected_tx = ""
    transcript_id = ""
    groupedTx = []
    for line in reader:
        #skip header line if it starts with '#'
        if not line[0].startswith('#'):
            #split rows by columns
            [chrom,program,source_or_exon,start_coord,stop_coord,tx_quality,strandedness,thephase,extra_annotation] = line            
            transcript_id = line[8].split(";")[1].split(" ")[2]
            #skip non-novel transcript, but save it...
            if any("reference_id" in attr for attr in line[8].split(";")):
                rejectTx.append(line)
                rejected_tx = transcript_id
                previous_tx = transcript_id
                continue
            #remove co-dependent line
            if rejected_tx == transcript_id:
                rejectTx.append(line)
                previous_tx = transcript_id
                continue
final_txs = []

#re-write new file with remaining correct transcripts:
              
#with open(path_to_files + 'a_subset.gtf','r') as bamfile:
with open(path_to_files + 'SRR1045766.bam.GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.harmonized_utrs.gtf','r') as bamfile:
    reader = csv.reader(bamfile, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in reader:
        #skip header line if it starts with '#'
        if not line[0].startswith('#'):
            if not line in rejectTx:
                [chrom,program,source_or_exon,start_coord,stop_coord,tx_quality,strandedness,thephase,extra_annotation] = line  
                transcript_id = line[8].split(";")[1].split(" ")[2]
                if source_or_exon == "exon":
                    start_coord, stop_coord = int(start_coord), int(stop_coord)
                    exons_dict[(chrom, strandedness, transcript_id)].append((start_coord, stop_coord))
                    introns_dict = defaultdict(set)
                    for tx_info, exons in exons_dict.items():
                        if len(exons) > 1:
                            [
                                    chrom,
                                    strandedness,
                                    tx_id
                                    ] = tx_info
                            exons = sorted(exons)
                            i=0
                            while (i + 1) < len(exons):
                                introns = (exons[i][1] + 1, exons[i+1][0] - 1)
                                i = i + 1                                
                                introns_dict[(chrom, strandedness, tx_id)].add(introns)

with gzip.open(path_to_files + '.gz','r') as fin:        
    for line in fin:        
        print('got line', line)
        #with open(path_to_files + '','r') as bamfile:
        #print(introns_dict)

#write_gtf_from_array(final_txs,path_to_files+"Novel_transcripts_with_introns.gtf")




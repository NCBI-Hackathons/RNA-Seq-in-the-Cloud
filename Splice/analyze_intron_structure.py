#Script by: SpliceItUp group
#UNC - Chapel Hill 
#Hackathon: RNA-seq-in-the-cloud

#Running this script locally:   python .py

import csv
import numpy as np
import pandas as pd

#path_to_files = "/home/vkkodali/example-gtf-files/"
path_to_files = "/home/ijimenez/RNA-Seq-in-the-Cloud/Splice/"
names=[]
ages=[]
rejectTx = []
intronless_Tx = []
#with open(path_to_files + 'SRR1045766.bam.GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.harmonized_utrs.gtf','r') as bamfile:
with open(path_to_files + 'mini_subset.gtf','r') as bamfile:
    reader = csv.reader(bamfile, delimiter = '\t')
    # check the header
    finding_intron = 0
    exon_counts = 0
    previous_tx = ""
    rejected_tx = ""
    transcript_id = ""
    for line in reader:
        #skip header line if it starts with '#'
        if not line[0].startswith('#'):
            #split rows by columns
            [chrom,program,source_or_exon,start_coord,stop_coord,tx_quality,strandedness,thephase,extra_annotation] = line            
            ages.append(strandedness)
            transcript_id = line[8].split(";")[2].split(" ")[2]
            print(transcript_id)
            if rejected_tx == line[8].split(";")[2].split(" ")[2]:
                rejectTx.append(line)
                previous_tx = transcript_id
                continue
            #skip non-novel transcript, but save it...
            if any("reference_id" in attr for attr in line[8].split(";")):
                rejectTx.append(line)
                rejected_tx = transcript_id
                previous_tx = transcript_id
                continue       
            if source_or_exon=="exon" and previous_tx == transcript_id:
                exon_counts = exon_counts+1
                print(exon_counts)
            else:
                exon_counts = 0
            if exon_counts <= 1:
                intronless_Tx.append(line)
                previous_tx = transcript_id
                continue
            previous_tx = transcript_id
#print(rejectTx)
print(intronless_Tx)




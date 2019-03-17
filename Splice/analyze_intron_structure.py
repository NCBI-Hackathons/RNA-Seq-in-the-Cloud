#Script by: SpliceItUp group
#UNC - Chapel Hill 
#Hackathon: RNA-seq-in-the-cloud

#Running this script locally:   python3 analyze_intron_structure.py

import io
import csv
import numpy as np
import pandas as pd
import glob
import os
import sys
import gzip
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL

path_to_files = "/data/study_gtf/"
#dir_path = os.path.dirname(os.path.realpath(__file__))
# path_to_AlignDB_reference="/home/ijimenez/aligndb_bams/"
# reference_tsv_file = path_to_AlignDB_reference+"combined_references.introns.ucsc.tsv.gz"
path_to_AlignDB_reference="/home/ijimenez/aligndb_bams/"
reference_tsv_file = "combined_references.introns.ucsc.tsv.gz"
#path_to_AlignDB_reference="/home/ijimenez/aligndb_bams/test_dir/"
#reference_tsv_file ="c_small_forced.gz"
path_to_output = "/home/ijimenez/output_test"
#path_to_test_file = path_to_files
#filename = "SRP009266_A549.gtf"
file_chunk = sys.argv[1]

import os
if not os.path.exists(path_to_output):
    os.makedirs(path_to_output)
output_path = path_to_output+"/"
def write_tsv_from_array(list_txs,file_name):
    find_unique_values = set(list_txs)
    the_list_of_txs = list(find_unique_values)
    with open(file_name, 'w') as file:
        for item in the_list_of_txs:
            item = '\t'.join(item)
            file.write("%s\n" % item)
        print("Created file! Name: "+file_name)

def find_introns_in_gtf(input_file_name):
    rejectTx = []
    intronless_Tx = []
    with open(path_to_files + input_file_name,'r') as gtf_file:
        reader = csv.reader(gtf_file, delimiter = '\t')
        # check the header
        finding_intron = 0
        exon_counts = 0
        previous_tx = ""
        rejected_tx = ""
        transcript_id = ""
        #groupedTx = []
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
    #final_txs = []
    with open(path_to_files + input_file_name,'r') as bamfile:
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
        introns_dict = {}
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
                    introns = str((exons[i][1] + 1, exons[i+1][0] - 1))
                    i = i + 1
                    introns_dict[(chrom, strandedness, tx_id)]=introns_dict.get((chrom,strandedness,tx_id),'')+','+introns
    # print("Introns from pseudo-novel transcripts:")
    return(introns_dict)

def find_introns_in_reference_tsv(ref_file_name):
    with gzip.open(path_to_AlignDB_reference + ref_file_name,'rt') as filein:      
        reader = csv.reader(filein, delimiter = '\t')
        ref_introns_dict = defaultdict(list)
        for actual_line in reader:
            if not actual_line[0].startswith('#') or not actual_line[0].startswith('W'):
                #ref_chrom,ref_transcriptid_unparsed,ref_strand
                line_copy = actual_line
                ref_introns = line_copy[-1]
                ref_introns = str(ref_introns)
                ref_chrom = line_copy[0]
                line_copy = actual_line
                ref_transcriptid_unparsed = line_copy[1]
                line_copy = actual_line
                ref_strand = line_copy[2]
                line_copy = actual_line
                array_of_tuples = ref_introns.split(",")
                if ref_strand == "-":
                    array_of_tuples = list(reversed(array_of_tuples))
                the_introns_key = ""
                #'(97364,99125),(88421,97293),(83318,88294),(80029,83212),(78953,79936),(58990,78841)': {('chrUn_GL000219v1', 'AK055615.1', '-')}
                for tuple in array_of_tuples:
                    [start_intron_pos, end_intron_pos] = tuple.split("..")
                    ref_start_pos = int(str(start_intron_pos).strip("[").strip("("))
                    ref_end_pos = int(str(end_intron_pos).strip("]").strip(")"))
                    if the_introns_key == "":
                        the_introns_key=the_introns_key+"("+str(ref_start_pos)+", "+str(ref_end_pos)+")"
                    else:
                        the_introns_key=the_introns_key+","+ "("+str(ref_start_pos)+", "+str(ref_end_pos)+")"
                ref_introns_dict[(the_introns_key)].append((ref_chrom,ref_strand))
                # print("Introns from reference transcripts:")
    return(ref_introns_dict)

#produce_tsv_transcripts_file
def produce_tsv_transcripts_file(intron_dict,reference_intron_dict):
    resulting_transcripts = []
    for tx_info,introns  in intron_dict.items():
        #print("GTF_Target: "+introns.lstrip(','))
        [chrom, strandedness, tx_id] = tx_info
        for ref_introns,ref_info in reference_intron_dict.items():
            #print(ref_info)
            if len(ref_info) > 1:
                for i in range(len(ref_info)):
                    #print(ref_info[i])
                    (ref_chrom, ref_strandedness) = ref_info[i]
                    #print("GOT EM: "+str(ref_chrom)+str(ref_strandedness))
                    if str(chrom) == str(ref_chrom):
                        if str(strandedness) == str(ref_strandedness):
                            #print("Comparing against ref: "+str(ref_introns))
                            if introns.lstrip(',') in ref_introns:
                                resulting_transcripts.append(tx_info)
            else:
                [(ref_chrom, ref_strandedness)] = ref_info
                #print("GOT EM: "+str(ref_chrom)+str(ref_strandedness))
                if str(chrom) == str(ref_chrom):
                    if str(strandedness) == str(ref_strandedness):
                        #print("Comparing against ref: "+str(ref_introns))
                        if introns.lstrip(',') in ref_introns:
                            resulting_transcripts.append(tx_info)
    return(resulting_transcripts)

the_reference_dictionary = find_introns_in_reference_tsv(reference_tsv_file)
#print("The reference dictionary: ")
#print(the_reference_dictionary)
#full_path_to_file = path_to_test_file+filename
#print("Running on "+filename+"...")
#split_the_path = full_path_to_file.split("/")
#the_input_file = split_the_path[-1]
#split_the_input_name = str(the_input_file).split(".gtf")
#the_output_file_name =split_the_input_name[0]
#the_output_file = "novel_transcripts_from_IntronStructure_"+the_output_file_name+".tsv" 
#the_input_introns_dictionary = find_introns_in_gtf(the_input_file)
#print("")
#print("TARGET DICTIONARY: ")
#print(the_input_introns_dictionary)
#write_tsv_from_array(produce_tsv_transcripts_file(the_input_introns_dictionary,the_reference_dictionary),output_path+the_output_file)

for filename in glob.iglob(path_to_files+file_chunk+'*.gtf'):
    full_path_to_file = filename
    print("Running on "+filename+"...")
    split_the_path = full_path_to_file.split("/")
    the_input_file = split_the_path[-1]
    split_the_input_name = str(the_input_file).split(".gtf")
    the_output_file_name =split_the_input_name[0]
    print(the_input_file)
    if the_input_file == "SRP009266_A549.gtf":
        print("Skipped the processed SRP file..."
        next
    the_output_file = "novel_transcripts_from_IntronStructure_"+the_output_file_name+".tsv" 
    the_input_introns_dictionary = find_introns_in_gtf(the_input_file)
    write_tsv_from_array(
        produce_tsv_transcripts_file(the_input_introns_dictionary,the_reference_dictionary),
        output_path+the_output_file
        )

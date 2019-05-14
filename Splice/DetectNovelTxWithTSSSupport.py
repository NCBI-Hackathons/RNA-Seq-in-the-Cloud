#!/usr/bin/env python3

import csv
import sys
import glob
import time
import os
import argparse
from ParseGTFFunctions import *

parser = argparse.ArgumentParser(description = """This script parses input gtf
    file along with CAGE clusters data in gff3 format to produce a gtf file with
    only novel transcripts supported by a CAGE cluster """)
parser.add_argument('-c', '--cage_clusters', help="cage clusters in gff3 \
    format; required")
parser.add_argument('-a', '--annotfile', help="annotation gtf; required")
args = parser.parse_args()

if not all([args.cage_clusters, args.annotfile]):
    print("Need annotation file and CAGE clusters file", file = sys.stderr)
    sys.exit()
else:
    annot_gtf = args.annotfile
    cage_gff3 = args.cage_clusters

start_time = time.time()

print("{0:4.2f} processing cage data..." .format(time.time() - start_time))
cage_clusters, cage_dict = process_cage_data(cage_gff3)
print('Total number of CAGE clusters: {}' .format(count_feats(cage_clusters)))

print("{0:4.2f} processing annot data..." .format(time.time() - start_time))
annot_tx = collect_novel_tx_from_exp_gtf(annot_gtf)
annot_first_exons, _ = collect_first_exons(annot_tx, min_exon_ct = 2)
print('Total number of first exons in the annotation: {}' .format(count_feats(annot_first_exons)))

gtf_files = [f for f in glob.glob('*.gtf') if f[:2] in ['SR','DR','ER']]

for exp_gtf in gtf_files:
    print("processing {}..." .format(exp_gtf))
    sra_acc = exp_gtf.split('.')[0]
    os.makedirs("./NovelTxWithTSSSupport/", exist_ok=True)
    cage_tx_gtf = './NovelTxWithTSSSupport/' + sra_acc + '.with_cage_support.gtf'

    print("{0:4.2f} collecting first exons of novel tx..." .format(time.time() - start_time))
    rnaseq_novel_tx = collect_novel_tx_from_exp_gtf(exp_gtf)
    rnaseq_first_exons, rnaseq_first_exon_tx_dict = collect_first_exons(rnaseq_novel_tx, min_exon_ct = 3)
    rnaseq_first_exon_ct = count_feats(rnaseq_first_exons)

    print("{0:4.2f} filter rnaseq exons sharing known donor sites..." .format(time.time() - start_time))
    rnaseq_first_exons = filter_first_exons_with_known_donors(rnaseq_first_exons, annot_first_exons)
    novel_first_exons_ct = count_feats(rnaseq_first_exons)

    print("{0:4.2f} mapping cage data..." .format(time.time() - start_time))
    first_exons_with_cage_support = find_tx_start_overlapping_cage(rnaseq_first_exons, cage_clusters)
    exons_with_cage_supp_ct = count_feats(first_exons_with_cage_support)

    print("{0:4.2f} writing output data..." .format(time.time() - start_time))
    write_filtered_gtf(exp_gtf, cage_tx_gtf, first_exons_with_cage_support, rnaseq_first_exon_tx_dict)

    print("Total number of first exons in the experimental GTF: {}" .format(rnaseq_first_exon_ct))
    print("Number of novel first exons in the experimental GTF: {}" .format(novel_first_exons_ct))
    print("Number of novel first exons with CAGE support: {}" .format(exons_with_cage_supp_ct))

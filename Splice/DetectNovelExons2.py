#!/usr/bin/env python3

import csv
import sys
import glob
import time
import os
import argparse
from ParseGTFFunctions import *

parser = argparse.ArgumentParser(description = """This script parses input gtf
    file along with annotation in gtf format and exons from transcript
    alignments in tsv format to produce a bed file with novel exons """)
parser.add_argument('-a', '--annotfile', help="annotation gtf; required")
parser.add_argument('-t', '--txexons', help="file with exons from transcript \
    alignments ; required")
args = parser.parse_args()

if not all([args.annotfile, args.txexons]):
    print("Need annotation file and tx exons file", file = sys.stderr)
    sys.exit()
else:
    annot_gtf = args.annotfile
    txalign_tsv = args.txexons

start_time = time.time()

print("{0:4.2f} collecting annot exons..." .format(time.time() - start_time))
annotated_exons = collect_exons_from_annotation(annot_gtf)

print("{0:4.2f} collecting tx alignment exons..." .format(time.time() - start_time))
aligndb_exons = collect_exons_from_tsv(txalign_tsv)

print("{0:4.2f} making a list of known exon starts and ends..." .format(time.time() - start_time))
known_exon_starts = defaultdict(set)
known_exon_ends = defaultdict(set)
process_exons_ends(annotated_exons, known_exon_starts, known_exon_ends)
# process_exons_ends(aligndb_exons, known_exon_starts, known_exon_ends)

gtf_files = [f for f in glob.glob('*.gtf') if f[:2] in ['SR','DR','ER']]

for exp_gtf in gtf_files:
    print("processing {}..." .format(exp_gtf))
    sra_acc = exp_gtf.split('.')[0]
    os.makedirs("./NovelExons/", exist_ok=True)
    unrep_exons_out = './NovelExons/' + sra_acc + '.unrep_exons.bed'
    novel_exons_out = './NovelExons/' + sra_acc + '.novel_exons.bed'

    print("{0:4.2f} collecting rnaseq_exons..." .format(time.time() - start_time))
    rnaseq_tx = collect_novel_tx_from_exp_gtf(exp_gtf)
    rnaseq_exons = collect_internal_exons(rnaseq_tx)
    total_rnaseq_exons_ct = count_feats(rnaseq_exons)

    print("{0:4.2f} removing annotated exons from rnaseq exons..." .format(time.time() - start_time))
    rnaseq_exons = collect_unique_feats(rnaseq_exons, annotated_exons)
    unannotated_exons_ct = count_feats(rnaseq_exons)

    print("{0:4.2f} removing retained introns..." .format(time.time() - start_time))
    rnaseq_exons = remove_retained_introns(rnaseq_exons, known_exon_starts, known_exon_ends)
    retained_introns_ct = unannotated_exons_ct - count_feats(rnaseq_exons)

    print("{0:4.2f} collecting novel exons..." .format(time.time() - start_time))
    novel_exons = collect_unique_feats(rnaseq_exons, aligndb_exons)
    novel_exons_ct = count_feats(novel_exons)

    print("{0:4.2f} collecting unrepresented exons..." .format(time.time() - start_time))
    unrep_exons = collect_common_feats(rnaseq_exons, aligndb_exons)
    unrep_exons_ct = count_feats(unrep_exons)

    print("writing files...")
    write_exons_to_file(novel_exons_out, novel_exons)
    write_exons_to_file(unrep_exons_out, unrep_exons)

    print("total number of internal exons: {}" .format(total_rnaseq_exons_ct))
    print("total number of exons after annot filtering: {}" .format(unannotated_exons_ct))
    print("total number of exons that retain introns: {}" .format(retained_introns_ct))
    print("total number of exons after aligndb filtering: {}" .format(novel_exons_ct))
    print("total number of unrepresented exons: {}" .format(unrep_exons_ct))

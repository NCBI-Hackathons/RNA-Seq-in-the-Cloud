#!/usr/bin/env python3

import csv
import sys
import glob
import time
import os
import argparse
from collections import defaultdict
from ParseGTFFunctions import *

parser = argparse.ArgumentParser(description = """This script parses input gtf
    file along with splice structures from transcript alignments in tsv format
    to produce a gtf file with unrepresented and novel transcripts """)
parser.add_argument('-a', '--annotfile', help="annotation gtf; required")
parser.add_argument('-t', '--txaligns', help="file with splice structures from \
    transcript alignments ; required")
args = parser.parse_args()

if not all([args.annotfile, args.txaligns]):
    print("Need annotation file and tx aligns file", file = sys.stderr)
    sys.exit()
else:
    annot_gtf = args.annotfile
    txalign_tsv = args.txaligns

start_time = time.time()

print("{0:4.2f} collecting splice structures of annot tx..." .format(time.time() - start_time), end='')
annot_tx = collect_novel_tx_from_exp_gtf(annot_gtf)
annot_tx = filter_tx_on_int_ct(annot_tx, min_int_count = 2)
annot_int_structures, _ = construct_splice_structures(annot_tx)
print('Done.')

print("{0:4.2f} collecting splice structures of alignments tx..." .format(time.time() - start_time), end='')
align_int_structures = process_tsv_splice_structures(txalign_tsv, min_int_count = 2)
print('Done.')

print('Input summary: ')
print('Total number of unique splice structures (with at least 2 introns) in the annotation: {}' .format(count_feats(annot_int_structures)))
print('Total number of unique splice structures (with at least 2 introns) found in tx alignments: {}\n' .format(count_feats(align_int_structures)))

gtf_files = [f for f in glob.glob('*.gtf') if f[:2] in ['SR','DR','ER']]

for exp_gtf in gtf_files:
    print("Summary for {}:" .format(exp_gtf))
    sra_acc = exp_gtf.split('.')[0]
    os.makedirs("./NovelTxWithIntronSupport/", exist_ok=True)
    unrep_int_file = './NovelTxWithIntronSupport/' + sra_acc + '.unrep_tx.gtf'
    novel_int_file = './NovelTxWithIntronSupport/' + sra_acc + '.novel_tx.gtf'

    print("{0:4.2f} collecting splice structures of novel tx..." .format(time.time() - start_time), end='')
    rnaseq_novel_tx = collect_novel_tx_from_exp_gtf(exp_gtf)
    rnaseq_novel_tx = filter_tx_on_int_ct(rnaseq_novel_tx, min_int_count = 2)
    rnaseq_int_structures, rnaseq_splice_struct_dict = construct_splice_structures(rnaseq_novel_tx)
    print('Done.')
    total_rnaseq_int_struct_ct = count_feats(rnaseq_int_structures)

    print("{0:4.2f} filtering annotated splice structures..." .format(time.time() - start_time), end='')
    rnaseq_int_structures = collect_unique_feats(rnaseq_int_structures, annot_int_structures)
    print('Done.')
    rnaseq_int_structs_not_in_annot = count_feats(rnaseq_int_structures)

    print("{0:4.2f} computing unrepresented splice structures..." .format(time.time() - start_time), end='')
    unrep_int_structures = collect_common_feats(rnaseq_int_structures, align_int_structures)
    print('Done.')
    unrep_int_struct_ct = count_feats(unrep_int_structures)

    print("{0:4.2f} writing unrepresented tx to gtf file..." .format(time.time() - start_time), end='')
    write_filtered_gtf(exp_gtf, unrep_int_file, unrep_int_structures, rnaseq_splice_struct_dict)
    print('Done.')

    print("{0:4.2f} novel splice structures..." .format(time.time() - start_time), end='')
    novel_int_structures = collect_unique_feats(rnaseq_int_structures, align_int_structures)
    print('Done.')
    novel_int_struct_ct = count_feats(novel_int_structures)

    print("{0:4.2f} writing novel tx to gtf file..." .format(time.time() - start_time), end='')
    write_filtered_gtf(exp_gtf, novel_int_file, novel_int_structures, rnaseq_splice_struct_dict)
    print('Done.\n')

    print('Total number of unique splice sturctures in {} file: {}' .format(exp_gtf, total_rnaseq_int_struct_ct))
    print('After removing splice structures that are already annotated, we have: {}' .format(rnaseq_int_structs_not_in_annot))
    print('Splice structures that are unrepresented: {}' .format(unrep_int_struct_ct))
    print('Novel splice structures: {}\n' .format(novel_int_struct_ct))

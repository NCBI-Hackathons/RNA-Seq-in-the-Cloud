#!/usr/bin/env python3

import gzip
import csv
import sys
import os
import re

if len(sys.argv) !=2 or sys.argv[1] in ['-h', '--h', '-help', '--help']:
    print("Usage: cat <input_file> | import_tx_alignment_exons.py <output_file>")
    sys.exit()
    
infile = sys.stdin
outfile = sys.argv[1]

tblin = csv.reader(infile, delimiter = '\t')
known_exon_set = set()
for line in tblin:
    if line[0].startswith('chr'):
        [chrom, tx, strand, exons] = line
        exons = exons.split(',')
        for exon in exons:
            exon = re.sub('[\[\]()]', '', exon)
            exon = tuple(exon.split('..'))
            known_exon_set.add((chrom, strand) + exon)
infile.close()

with gzip.open(outfile, 'wt') as f:
    tbl = csv.writer(f, delimiter = '\t', lineterminator = os.linesep)
    for known_exon in sorted(known_exon_set):
        tbl.writerow(known_exon)

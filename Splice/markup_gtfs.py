#!/usr/bin/env python3

import csv
import sys
import os
import re

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'rt') as fi, open(outfile, 'wt') as fo:
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(fo,
            delimiter = '\t',
            lineterminator = os.linesep,
            quotechar = '\xb6' )
    for line in tblin:
        if not line[0].startswith('#'):
            if 'reference_id' not in line[8]:
                tx_id = re.search('transcript_id "STRG.(\d+).\d+"', line[8]).group(1)
                if int(tx_id) % 3 == 0:
                    line[8] = line[8] + 'novel "yes";'
        tblout.writerow(line)

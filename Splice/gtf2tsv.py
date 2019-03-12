#!/usr/bin/env python

import csv
import gzip
import re
import sys
import os
from collections import defaultdict, OrderedDict
import argparse

# input_file = '/home/vamsi/Downloads/sample.gtf'
# output_file = '/home/vamsi/Downloads/sample.tsv'

parser = argparse.ArgumentParser(description ="""This script parses input
                gtf file and returns a tab-delimited file with one column
                for each attribute in column 9 """)
parser.add_argument('-i', '--infile', help="input file; required")
parser.add_argument('-o', '--outfile', help="output file; required")
args = parser.parse_args()

if not all([args.infile, args.outfile]):
    print("You must provide input and output files", file = sys.stderr)
    sys.exit()
else:
    input_file = args.infile
    output_file = args.outfile

def create_attrib_dict(input_file):
    attrib_list = set()
    with open(input_file, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            if not line[0].startswith('#'):
                attribs = process_attribs(line[8])
                for attrib in attribs.keys():
                    attrib_list.add(attrib)
    return attrib_list

def process_attribs(attribs):
    new_attribs = {}
    attribs = filter(None, attribs.rstrip(' ').split(';'))
    for attrib in attribs:
        attrib = attrib.lstrip(' ').split(' ')
        if attrib[0] == 'db_xref':
            attrib = OrderedDict.fromkeys(attrib[1].strip('"').split(':'))
            attrib = list(attrib.keys())
        new_attribs[attrib[0]] = attrib[1].strip('"')
    return new_attribs

all_attribs = create_attrib_dict(input_file)
all_attribs = list(all_attribs)

with open(input_file, 'rt') as fi, open(output_file, 'wt') as fo:
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(fo, delimiter = '\t', lineterminator = os.linesep)
    tblout.writerow(['chrom', 'source', 'feat_type', 'start', 'stop',
                    'score', 'strand', 'phase'] + all_attribs)
    for line in tblin:
        if not line[0].startswith('#'):
            attribs = process_attribs(line[8])
            attrib_vals = [attribs.get(i, 'NULL') for i in all_attribs]
            tblout.writerow(line[:8] + attrib_vals)

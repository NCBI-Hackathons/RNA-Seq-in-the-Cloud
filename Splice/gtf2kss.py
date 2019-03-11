#!/usr/bin/env python

import csv
import os
import sys
import gzip
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description ="""This script parses input
                gtf file and generates a table of splice sites to be used
                with `--known-splice-sites` parameter of hisat2""")
parser.add_argument('-i', '--infile', help="input file; default STDIN")
parser.add_argument('-o', '--outfile', help="output file; default STDOUT")
parser.add_argument('-f', '--format', help="input file format; default is gtf")
parser.add_argument('-a', '--assembly', help="assembly_report.txt file")
parser.add_argument('-e', '--exons', help="filename for exons file")
parser.add_argument('-z', '--gzip', help="input is gzip compressed", \
                    action = 'store_true')
args = parser.parse_args()

if args.infile:
    if args.gzip:
        input_file = gzip.open(args.infile, 'rt')
    else:
        input_file = open(args.infile, 'r')
else:
    input_file = sys.stdin

if args.outfile:
    splice_sites_file = open(args.outfile, 'w')
else:
    splice_sites_file = sys.stdout

if args.format:
    input_format = args.format
    if not input_format in ['gff3', 'gtf']:
        print("Input file format can only be `gff3` or `gtf`")
        sys.exit()
else:
    input_format = 'gtf'

if args.assembly:
    assembly_report_file = args.assembly
    with open(assembly_report_file, 'r') as f:
        line = f.readline()
        if not line.startswith('# Assembly name:'):
            print("Assembly report file is not correctly formatted. Please "
                  "download the assembly_report.txt file from NCBI asssembly "
                  "for the assembly of choice")
            sys.exit()

def create_mapping_dict(assembly_report_file):
    """ returns a dictonary of seq-ids mapped to UCSC names;
    requires an assembly_report.txt file """
    mapping_dict = {}
    with open(assembly_report_file, 'r') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            if not line[0].startswith('#'):
                [
                    ens_name,
                    _,
                    _,
                    _,
                    gb_acc,
                    _,
                    rs_acc,
                    _,
                    _,
                    ucsc_name
                ] = line
                mapping_dict[ens_name] = ucsc_name
                mapping_dict[gb_acc] = ucsc_name
                mapping_dict[rs_acc] = ucsc_name
                mapping_dict['CHR_' + ens_name] = ucsc_name
                mapping_dict[ucsc_name] = ucsc_name
    return mapping_dict

def get_gtf_exons(input_file):
    """ returns exons for each transcript in a gtf file;
    requires gtf file as input """
    def process_attribs(attribs):
        new_attribs = {}
        attribs = filter(None, attribs.rstrip(' ').split(';'))
        for attrib in attribs:
            attrib = attrib.lstrip(' ').split(' ')
            new_attribs[attrib[0]] = attrib[1].strip('"')
        return new_attribs

    tbl = csv.reader(input_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom,
                feat_source,
                feat_type,
                start,
                stop,
                score,
                strand,
                phase,
                attribs
            ] = line
            start, stop = int(start), int(stop)
            if feat_type == "exon" and stop > start:
                new_attribs = process_attribs(attribs)
                tx = new_attribs['transcript_id']
                exons_dict[(chrom, strand, tx)].append([start, stop])
    input_file.close()
    return exons_dict

def get_gff_exons(input_file):
    """ returns exons for each transcript in a gff3 file;
    requires gff3 file as input """
    def process_attribs(attribs):
        new_attribs = {}
        attribs = list(filter(None, attribs.split(';'))) ## removes empty strings, needed because some gff3 lines have ";;"
        for attrib in attribs:
            k, v = attrib.split('=')
            if k == 'Dbxref':
                xrefs = v.split(',')
                for xref in xrefs:
                    terms = xref.split(':')
                    new_attribs[terms[-2]] = terms[-1]
            else:
                new_attribs[k] = v
        return new_attribs

    tbl = csv.reader(input_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom,
                feat_source,
                feat_type,
                start,
                stop,
                score,
                strand,
                phase,
                attribs
            ] = line
            start, stop = int(start), int(stop)
            if feat_type == 'exon' and stop > start:
                new_attribs = process_attribs(attribs)
                if 'transcript_id' in attribs:
                    tx = new_attribs['transcript_id']
                else:
                    tx = new_attribs['Parent']
                exons_dict[(chrom, strand, tx)].append([start, stop])
    input_file.close()
    return exons_dict

def get_splice_sites_and_exons(exons_dict):
    """ generates a list of unique splice sites;
    requires a dictionary of exons """
    splice_sites = set()
    unique_exons = defaultdict(set)
    for rs_info, exons in exons_dict.items():
        [
            chrom,
            strand,
            rs
        ] = rs_info
        exons.sort()
        tx_exons = [exons[0]]
        i = 0
        while (i + 1) < len(exons):
            if exons[i+1][0] - tx_exons[-1][1] > 5:
                tx_exons.append(exons[i+1])
            else:
                tx_exons[-1][1] = exons[i+1][1]
            i = i + 1
        if len(tx_exons) > 1:
            for i in range(1, len(tx_exons)):
                splice_sites.add((chrom,
                                tx_exons[i-1][1] - 1,
                                tx_exons[i][0] - 1,
                                strand)) ## 0-based coordinates
        for exons in tx_exons:
            unique_exons[(chrom, strand)].add((exons[0], exons[1]))

    for (chrom, strand), exons in unique_exons.items():
        exons = sorted(list(exons))
        merged_exons = [list(exons[0])]
        for exon in exons:
            exon = list(exon)
            if exon[0] == merged_exons[-1][0]:
                merged_exons[-1] = exon
            elif exon[1] == merged_exons[-1][1]:
                continue
            elif exon[0] < merged_exons[-1][1]:
                merged_exons[-1][1] = exon[1]
            else:
                merged_exons.append(exon)
        unique_exons[(chrom, strand)] = merged_exons

    merged_exons = set()
    for (chrom, strand), exons in unique_exons.items():
        for exon in exons:
            merged_exons.add((chrom, exon[0]-1, exon[1]-1, strand))

    return splice_sites, merged_exons

def remap_seqids(feat_set, mapping_dict):
    remapped_feats = set()
    for feat in feat_set:
        [chrom, start, stop, strand] = feat
        if mapping_dict[chrom] != 'na':
            chrom = mapping_dict[chrom]
        remapped_feats.add((chrom, start, stop, strand))
    return remapped_feats


def tabulate_exons(unique_exons, exons_file):
    with open(exons_file, 'w') as f:
        tbl = csv.writer(f,
                        delimiter = '\t',
                        lineterminator = os.linesep)
        tbl.writerow(['#chromosome',
                      'exon_start',
                      'exon_end',
                      'strand'])
        tbl.writerows(sorted(list(unique_exons)))

def tabulate_splice_sites(splice_sites, splice_sites_file):
    tbl = csv.writer(splice_sites_file,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#chromosome',
                  'splice_donor',
                  'splice_acceptor',
                  'strand'])
    tbl.writerows(sorted(list(splice_sites)))
    splice_sites_file.close()

if input_format == 'gff3':
    exons_dict = get_gff_exons(input_file)
elif input_format == 'gtf':
    exons_dict = get_gtf_exons(input_file)
splice_sites, unique_exons = get_splice_sites_and_exons(exons_dict)

if args.assembly:
    mapping_dict = create_mapping_dict(assembly_report_file)
    mapped_splice_sites = remap_seqids(splice_sites, mapping_dict)
    tabulate_splice_sites(mapped_splice_sites, splice_sites_file)
else:
    tabulate_splice_sites(splice_sites, splice_sites_file)

if args.exons:
    exons_file = args.exons
    if args.assembly:
        mapped_unique_exons = remap_seqids(unique_exons, mapping_dict)
        tabulate_exons(mapped_unique_exons, exons_file)
    else:
        tabulate_exons(unique_exons, exons_file)

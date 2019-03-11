#!/usr/bin/env python

import csv
import gzip
import re
import sys
import os
from collections import defaultdict, OrderedDict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description ="""This script parses input
                gff/gtf file and harmonizes the equivalent terminal exons
                for each gene """)
parser.add_argument('-i', '--infile', help="input file; required")
parser.add_argument('-o', '--outfile', help="output file; required")
parser.add_argument('-f', '--format', help="input file format `gff3` or `gtf`; \
    default is gtf")
parser.add_argument('--fancy_refseq_stuff', help="do fancy refseq stuff to \
    extend only up to the longest known RefSeq on the 5' end and to the \
    longest 'any' RefSeq on the 3' end", action = 'store_true')
parser.add_argument('-z', '--gzip', help="input and ouput files are gzip \
    compressed", action = 'store_true')
args = parser.parse_args()

if not all([args.infile, args.outfile]):
    print("You must provide input and output files", file = sys.stderr)
    sys.exit()
else:
    input_file = args.infile
    output_file = args.outfile

if args.format:
    input_format = args.format
    if not input_format in ['gff3', 'gtf']:
        print("Input file format can only be `gff3` or `gtf`", file = sys.stderr)
        sys.exit()
else:
    input_format = 'gtf'

if args.fancy_refseq_stuff:
    fancy_refseq_stuff = True
else:
    fancy_refseq_stuff = False

def get_gff3_exons(input_file):
    """ returns the first exon for each transcript in a gff3 file;
    requires gff3 file as input """
    if args.gzip:
        f = gzip.open(input_file, 'rt')
    else:
        f = open(input_file, 'rt')
    def process_attribs(attribs):
        new_attribs = {}
        attribs = list(filter(None, attribs.split(';')))
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

    tbl = csv.reader(f, delimiter = '\t')
    parent_id_dict = {}
    first_exons_dict = {}
    last_exons_dict = {}
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom, feat_source, feat_type,
                start, stop, score,
                strand, phase, attribs
            ] = line
            start, stop = int(start), int(stop)
            if feat_type == 'exon' and stop > start:
                new_attribs = process_attribs(attribs)
                gene_id = new_attribs['GeneID']
                if 'transcript_id' in attribs:
                    tx = new_attribs['transcript_id']
                else:
                    tx = new_attribs['Parent']
                gene_info = (chrom, strand, gene_id)
                parent_id_dict[new_attribs['Parent']] = (
                    chrom, strand, gene_id, tx)
                if gene_info not in first_exons_dict:
                    first_exons_dict[gene_info] = {tx : [start, stop]}
                    last_exons_dict[gene_info] = {tx : [start, stop]}
                elif tx not in first_exons_dict[gene_info]:
                    first_exons_dict[gene_info][tx] = [start, stop]
                    last_exons_dict[gene_info][tx] = [start, stop]
                elif strand == '+':
                    if start < first_exons_dict[gene_info][tx][0]:
                        first_exons_dict[gene_info][tx] = [start, stop]
                    if start > last_exons_dict[gene_info][tx][0]:
                        last_exons_dict[gene_info][tx] = [start, stop]
                elif strand == '-':
                    if start > first_exons_dict[gene_info][tx][0]:
                        first_exons_dict[gene_info][tx] = [start, stop]
                    if start < last_exons_dict[gene_info][tx][0]:
                        last_exons_dict[gene_info][tx] = [start, stop]
    f.close()
    return first_exons_dict, last_exons_dict, parent_id_dict

def get_gtf_exons(input_file):
    """ returns first exon for each transcript in a gtf file;
    requires gtf file as input
    process_attribs notes (mostly for my future self):
    1. there can be more attribs like db_xref that are present >1 time; this
    script ends up keeping only the one that appears in the end
    2. 'description' attrib will include only the first word
    3. used OrderedDict b/c list->set->list conversion may lose order and we
    will end up with GeneID:123 in one case and 123:GeneID in another
    """
    if args.gzip:
        f = gzip.open(input_file, 'rt')
    else:
        f = open(input_file, 'rt')
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

    tbl = csv.reader(f, delimiter = '\t')
    first_exons_dict = {}
    last_exons_dict = {}
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom, feat_source, feat_type,
                start, stop, score,
                strand, phase, attribs
            ] = line
            start, stop = int(start), int(stop)
            if feat_type == 'exon' and stop > start:
                new_attribs = process_attribs(attribs)
                gene_id = new_attribs.get('GeneID', new_attribs['gene_id'])
                tx = new_attribs['transcript_id']
                gene_info = (chrom, strand, gene_id)
                if gene_info not in first_exons_dict:
                    first_exons_dict[gene_info] = {tx : [start, stop]}
                    last_exons_dict[gene_info] = {tx : [start, stop]}
                elif tx not in first_exons_dict[gene_info]:
                    first_exons_dict[gene_info][tx] = [start, stop]
                    last_exons_dict[gene_info][tx] = [start, stop]
                elif strand == '+':
                    if start < first_exons_dict[gene_info][tx][0]:
                        first_exons_dict[gene_info][tx] = [start, stop]
                    if start > last_exons_dict[gene_info][tx][0]:
                        last_exons_dict[gene_info][tx] = [start, stop]
                elif strand == '-':
                    if start > first_exons_dict[gene_info][tx][0]:
                        first_exons_dict[gene_info][tx] = [start, stop]
                    if start < last_exons_dict[gene_info][tx][0]:
                        last_exons_dict[gene_info][tx] = [start, stop]
    f.close()
    return first_exons_dict, last_exons_dict

def create_start_extensions_list(first_exons_dict, fancy_refseq_stuff):
    start_extensions = {}
    for gene_info, all_txs in first_exons_dict.items():
        chrom, strand, gene_id = gene_info
        longest_first_exons = {}
        if fancy_refseq_stuff:
            if strand == '+':
                for tx, first_exon in sorted(all_txs.items()):
                    i = longest_first_exons.get(first_exon[1], None)
                    if tx.startswith('N'):
                        if i is None or first_exon[0] < i :
                            longest_first_exons[first_exon[1]] = first_exon[0]
                    elif first_exon[1] not in longest_first_exons:
                        if i is None or first_exon[0] < i :
                            longest_first_exons[first_exon[1]] = first_exon[0]
            elif strand == '-':
                for tx, first_exon in sorted(all_txs.items()):
                    i = longest_first_exons.get(first_exon[0], None)
                    if tx.startswith('N'):
                        if i is None or first_exon[1] > i:
                            longest_first_exons[first_exon[0]] = first_exon[1]
                    elif first_exon[0] not in longest_first_exons:
                        if i is None or first_exon[1] > i:
                            longest_first_exons[first_exon[0]] = first_exon[1]
        else:
            first_exons_set = set()
            for tx, first_exon in sorted(all_txs.items()):
                first_exons_set.add(tuple(first_exon))
            first_exons_set = sorted(list(first_exons_set))
            if strand == '+':
                for first_exon in first_exons_set:
                    if first_exon[1] not in longest_first_exons:
                        longest_first_exons[first_exon[1]] = first_exon[0]
            elif strand == '-':
                for first_exon in first_exons_set:
                    longest_first_exons[first_exon[0]] = first_exon[1]

        for tx, first_exon in all_txs.items():
            tx_info = (chrom, strand, gene_id, tx)
            if strand == '+':
                if longest_first_exons[first_exon[1]] == first_exon[0]:
                    continue
                else:
                    start_extensions[tx_info] = [
                         first_exon,
                         [longest_first_exons[first_exon[1]],
                         first_exon[1]]]
            elif strand == '-':
                if longest_first_exons[first_exon[0]] == first_exon[1]:
                    continue
                else:
                    start_extensions[tx_info] = [
                         first_exon,
                         [first_exon[0],
                         longest_first_exons[first_exon[0]]]]
    return start_extensions

def create_end_extensions_list(last_exons_dict):
    end_extensions = {}
    for gene_info, all_txs in last_exons_dict.items():
        chrom, strand, gene_id = gene_info

        last_exons_set = set()
        for tx, last_exon in all_txs.items():
            last_exons_set.add(tuple(last_exon))

        longest_last_exons = {}
        last_exons_set = sorted(list(last_exons_set))
        if strand == '+':
            for last_exon in last_exons_set:
                longest_last_exons[last_exon[0]] = last_exon[1]
        elif strand == '-':
            for last_exon in last_exons_set:
                if last_exon[1] not in longest_last_exons:
                    longest_last_exons[last_exon[1]] = last_exon[0]

        for tx, last_exon in all_txs.items():
            tx_info = (chrom, strand, gene_id, tx)
            if strand == '+':
                if longest_last_exons[last_exon[0]] == last_exon[1]:
                    continue
                else:
                    end_extensions[tx_info] = [
                         last_exon,
                         [last_exon[0],
                         longest_last_exons[last_exon[0]]]]
            elif strand == '-':
                if longest_last_exons[last_exon[1]] == last_exon[0]:
                    continue
                else:
                    end_extensions[tx_info] = [
                         last_exon,
                         [longest_last_exons[last_exon[1]],
                         last_exon[1]]]
    return end_extensions

def write_updated_gff3_file(input_file, output_file, start_extensions, end_extensions, parent_id_dict):
    if args.gzip:
        fi = gzip.open(input_file, 'rt')
        fo = gzip.open(output_file, 'wt')
    else:
        fi = open(input_file, 'rt')
        fo = open(output_file, 'rt')

    def process_attribs(attribs):
        new_attribs = {}
        split_attribs = list(filter(None, attribs.split(';')))
        for attrib in split_attribs:
            k, v = attrib.split('=')
            if k == 'Dbxref':
                xrefs = v.split(',')
                for xref in xrefs:
                    terms = xref.split(':')
                    new_attribs[terms[-2]] = terms[-1]
            else:
                new_attribs[k] = v
        return new_attribs

    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
        fo,
        delimiter = '\t',
        lineterminator = os.linesep)
    for line in tblin:
        if not line[0].startswith('#'):
            [
                chrom, feat_source, feat_type,
                start, stop, score,
                strand, phase, attribs
            ] = line
            start, stop = int(start), int(stop)
            new_attribs = process_attribs(attribs)
            if feat_type == 'exon' and stop > start:
                gene_id = new_attribs['GeneID']
                if 'transcript_id' in attribs:
                    tx = new_attribs['transcript_id']
                else:
                    tx = new_attribs['Parent']
                gene_info = (chrom, strand, gene_id, tx)
                if (gene_info in start_extensions
                   and start_extensions[gene_info][0] == [start, stop]):
                    line[3], line[4] = start_extensions[gene_info][1]
                elif (gene_info in end_extensions
                     and end_extensions[gene_info][0] == [start, stop]):
                    line[3], line[4] = end_extensions[gene_info][1]
            else:
                gff3_id = new_attribs.get('ID', 'No_ID_present')
                if (gff3_id in parent_id_dict
                   and parent_id_dict[gff3_id] in start_extensions):
                    gene_info = parent_id_dict[gff3_id]
                    if (strand == '+'
                       and start_extensions[gene_info][0][0] == start):
                        line[3] = start_extensions[gene_info][1][0]
                        line[8] = line[8].rstrip() + ';tag=extended_start'
                    elif (strand == '-'
                         and start_extensions[gene_info][0][1] == stop):
                        line[4] = start_extensions[gene_info][1][1]
                        line[8] = line[8].rstrip() + ';tag=extended_start'
                if (gff3_id in parent_id_dict
                   and parent_id_dict[gff3_id] in end_extensions):
                    gene_info = parent_id_dict[gff3_id]
                    if (strand == '+'
                       and end_extensions[gene_info][0][1] == stop):
                        line[4] = end_extensions[gene_info][1][1]
                        if line[8].endswith('extended_start'):
                            line[8] = line[8].rstrip() + ',extended_end'
                        else:
                            line[8] = line[8].rstrip() + ';tag=extended_end'
                    elif (strand == '-'
                         and end_extensions[gene_info][0][0] == start):
                        line[3] = end_extensions[gene_info][1][0]
                        if line[8].endswith('extended_start'):
                            line[8] = line[8].rstrip() + ',extended_end'
                        else:
                            line[8] = line[8].rstrip() + ';tag=extended_end'
        tblout.writerow(line)
    fi.close()
    fo.close()

def write_updated_gtf_file(input_file, output_file, start_extensions, end_extensions):
    if args.gzip:
        fi = gzip.open(input_file, 'rt')
        fo = gzip.open(output_file, 'wt')
    else:
        fi = open(input_file, 'rt')
        fo = open(output_file, 'wt')
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

    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(fo,
            delimiter = '\t',
            lineterminator = os.linesep,
            quotechar = '\xb6' )
    for line in tblin:
        if not line[0].startswith('#'):
            [
                chrom, feat_source, feat_type,
                start, stop, score,
                strand, phase, attribs
            ] = line
            start, stop = int(start), int(stop)
            if feat_type == 'exon' and stop > start:
                new_attribs = process_attribs(attribs)
                gene_id = new_attribs.get('GeneID', new_attribs['gene_id'])
                tx = new_attribs['transcript_id']
                gene_info = (chrom, strand, gene_id, tx)
                if (gene_info in start_extensions
                   and start_extensions[gene_info][0] == [start, stop]):
                    line[3], line[4] = start_extensions[gene_info][1]
                elif (gene_info in end_extensions
                     and end_extensions[gene_info][0] == [start, stop]):
                    line[3], line[4] = end_extensions[gene_info][1]
        tblout.writerow(line)
    fi.close()
    fo.close()

if input_format == 'gff3':
    first_exons_dict, last_exons_dict, parent_id_dict = get_gff3_exons(input_file)
else:
    first_exons_dict, last_exons_dict = get_gtf_exons(input_file)

start_extensions = create_start_extensions_list(first_exons_dict, fancy_refseq_stuff)
end_extensions = create_end_extensions_list(last_exons_dict)

if input_format == 'gff3':
    write_updated_gff3_file(
        input_file,
        output_file,
        start_extensions,
        end_extensions,
        parent_id_dict )
else:
    write_updated_gtf_file(
        input_file,
        output_file,
        start_extensions,
        end_extensions )

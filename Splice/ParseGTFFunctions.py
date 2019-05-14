import sys
import os
import re
import gzip
import csv
from collections import defaultdict, OrderedDict

def collect_novel_tx_from_exp_gtf(gtffile):
    """
    parse `gtffile` produced by StringTie to extract all novel transcripts
    that have at least one intron
    """
    novel_transcripts = defaultdict(list)
    with open(gtffile, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            if not line[0].startswith('#'):
                [ chrom, source, feat_type,
                  start, stop, score,
                  strand, phase, attribs ] = line
                if ('reference_id' not in attribs
                    and feat_type == 'exon'):
                    attribs = process_attribs(attribs)
                    gene_id = attribs['gene_id']
                    tx_id = attribs['transcript_id']
                    tx_info = (chrom, strand, gene_id, tx_id)
                    exon = (int(start), int(stop))
                    novel_transcripts[tx_info].append(exon)
    return novel_transcripts

def collect_internal_exons(novel_transcripts):
    """
    parse the novel_transcripts dict to extract a unique set of exons
    """
    internal_exons = defaultdict(set)
    for tx_info, exons in novel_transcripts.items():
        chrom, strand, gene_id, tx_id = tx_info
        if len(exons) > 2:
            exons.sort()
            for exon in exons[1:-1]:
                internal_exons[(chrom, strand)].add(exon)
    return internal_exons

def process_attribs(attribs):
    new_attribs = {}
    attribs = filter(None, attribs.rstrip(' ').split(';'))
    for attrib in attribs:
        attrib = re.split(' |=', attrib.lstrip(' '))
        if attrib[0] == 'db_xref':
            attrib = OrderedDict.fromkeys(attrib[1].strip('"').split(':'))
            attrib = list(attrib.keys())
        new_attribs[attrib[0]] = attrib[1].strip('"')
    return new_attribs

def filter_tx_on_int_ct(tx_dict, min_int_count = 2):
    """
    Parse trnascript dict to extract tx with introns >= min_int_count
    """
    filtered_tx = {}
    for tx_info, exons in tx_dict.items():
        if len(exons) >= min_int_count:
            filtered_tx[tx_info] = exons
    return filtered_tx

def process_cage_data(cage_gff3):
    """
    parses the GFF3 file with CAGE data and extracts the CAGE cluster feats
    """
    cage_clusters = defaultdict(set)
    cage_dict = {}
    with open(cage_gff3, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            if not line[0].startswith('#'):
                [ chrom, source, feat,
                start, stop, score,
                strand, phase, attribs ] = line
                if feat == 'region':
                    attribs = process_attribs(attribs)
                    cluster = (int(start), int(stop))
                    cluster_info = (attribs['standard_name'],
                        attribs['num_tags_in_cluster'])
                    cage_clusters[(chrom, strand)].add(cluster)
                    cage_dict[(chrom, strand) + cluster] = cluster_info
    return cage_clusters, cage_dict

def collect_first_exons(tx_dict, min_exon_ct = 2):
    """
    parses a dictionary object with novel transcripts and returns their first
    exons and a dictionary mapping tx names to those first exons.
    """
    first_exons = defaultdict(set)
    first_exon_tx_dict = defaultdict(set)
    for tx_info, exons in tx_dict.items():
        if len(exons) >= min_exon_ct:
            chrom, strand, gene_id, tx_id = tx_info
            exons.sort()
            i = 0 if strand == '+' else -1
            first_exons[(chrom, strand)].add(exons[i])
            if (chrom, strand) not in first_exon_tx_dict:
                first_exon_tx_dict[(chrom, strand)] = defaultdict(list)
            first_exon_tx_dict[(chrom, strand)][exons[i]].append(tx_info[2:])
    return first_exons, first_exon_tx_dict

def collect_exons_from_annotation(gtffile):
    """
    `gtffile` needs to be an uncompressed gtf file
    """
    chr_whitelist = [
    'chr1', 'chr10', 'chr11', 'chr12', 'chr13',
    'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
    'chr19', 'chr2', 'chr20', 'chr21', 'chr22',
    'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
    'chr8', 'chr9', 'chrM', 'chrX', 'chrY' ]
    gtf_exons = defaultdict(set)
    with open(gtffile, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            if (not line[0].startswith('#')
                and line[2] == 'exon'
                and line[0] in chr_whitelist ):
                chrom, strand = line[0], line[6]
                exon = (int(line[3]), int(line[4]))
                gtf_exons[(chrom, strand)].add(exon)
    return gtf_exons

def collect_exons_from_tsv(exons_tsv):
    """
    `exons_tsv` file should be a gzip-compressed, tab-delimited file with the
    first four fields chrom, strand, exon start and exon end.
    All co-ordinates are 1-based
    """
    tsv_exons = defaultdict(set)
    with gzip.open(exons_tsv, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        for line in tbl:
            tsv_exons[(line[0], line[1])].add((int(line[2]), int(line[3])))
    return tsv_exons

def process_exons_ends(all_exons, exon_starts, exon_ends):
    """
    need to provide `exon_starts` and `exon_ends` objects that are defaultdicts
    """
    for (chrom, strand), exons in all_exons.items():
        for exon in exons:
            exon_starts[(chrom, strand)].add(exon[0])
            exon_ends[(chrom, strand)].add(exon[1])
    return exon_starts, exon_ends

def remove_retained_introns(all_exons, exon_starts, exon_ends):
    """
    provided a dict with exons (`all_exons`) along with known starts
    (`exon_starts`) and ends (`exon_ends`), remove exons that retain introns
    """
    filtered_exons = defaultdict(set)
    for (chrom, strand), exons in all_exons.items():
        for exon in exons:
            exon_start, exon_end = exon
            if (exon_start in exon_starts[(chrom, strand)]
                and exon_end in exon_ends[(chrom, strand)]):
                continue
            else:
                filtered_exons[(chrom, strand)].add(exon)
    return filtered_exons

def collect_unique_feats(all_feats, known_feats):
    """
    `all_feats` and `known_feats` are defaultdict(set) type objects with data in
    the following format:
    {('chr1','+'):{(1,10),(20,30),(40,50)}, ('chr2', '-'):{(11,19),(31,42)}}
    returns another object similar to them called `filtered_feats` that
    contains exons unique to `all_exons`
    """
    unique_feats = {}
    for (chrom, strand), feats in all_feats.items():
        feats = feats - known_feats[(chrom, strand)]
        unique_feats[(chrom, strand)] = feats
    return unique_feats

def collect_common_feats(all_feats, known_feats):
    """
    `all_feats` and `known_feats` are defaultdict(set) type objects with data in
    the following format:
    {('chr1','+'):{(1,10),(20,30),(40,50)}, ('chr2', '-'):{(11,19),(31,42)}}
    returns another object similar to them called `filtered_feats` with
    exons common to both.
    """
    common_feats = {}
    for (chrom, strand), feats in all_feats.items():
        feats = feats & known_feats[(chrom, strand)]
        common_feats[(chrom, strand)] = feats
    return common_feats

def count_feats(feat_dict):
    """
    `feat_dict` is a defaultdict(set) type object with data in the following
    format: {('chr1','+'):{(1,10),(20,30),(40,50)}, ('chr2', '-'):{(11,19),(31,42)}}
    """
    count = 0
    for info, feats in feat_dict.items():
        count = count + len(feats)
    return count

def write_exons_to_file(outfile, all_exons):
    """
    writes all exons present in the `all_exons` dict to `outfile`
    `all_exons` is a dict with data in the following format:
    {('chr1','+'):{(1,10),(20,30),(40,50)}, ('chr2', '-'):{(11,19),(31,42)}}
    """
    with open(outfile, 'wt') as f:
        tbl = csv.writer(f, delimiter = '\t')
        for (chrom, strand), exons in all_exons.items():
            for exon in exons:
                tbl.writerow([chrom, exon[0]-1, exon[1], "novel_exon", "1000", strand])

def find_tx_start_overlapping_cage(first_exons, cage_clusters, dist = 0):
    clusters_with_overlapping_tx_starts = defaultdict(set)
    for (chrom, strand), all_first_exons in first_exons.items():
        sorted_first_exons = list(sorted(all_first_exons))
        sorted_cage_clusters = list(sorted(cage_clusters[(chrom, strand)]))
        i = 0 if strand == '+' else 1
        for exon in sorted_first_exons:
            for cluster in sorted_cage_clusters:
                if (exon[i] >= cluster[0] - dist
                    and exon[i] <= cluster[1] + dist ):
                    clusters_with_overlapping_tx_starts[(chrom, strand)].add(exon)
                    del sorted_cage_clusters[:sorted_cage_clusters.index(cluster)]
                    break
    return clusters_with_overlapping_tx_starts

def filter_first_exons_with_known_donors(all_first_exons, known_first_exons):
    filtered_first_exons = defaultdict(set)
    for (chrom, strand), all_exons in all_first_exons.items():
        i = 1 if strand == '+' else 0
        known_donors = set(exon[i] for exon in known_first_exons[(chrom, strand)])
        for exon in all_exons:
            if exon[i] in known_donors:
                continue
            else:
                filtered_first_exons[(chrom, strand)].add(exon)
    return filtered_first_exons

def construct_splice_structures(tx_dict):
    """
    Parse `tx_dict` dictionary of transcripts to generate a set of unique
    splice structures for each chrom/strand and a dictionary of
    splice_structure : tx_info.
    """
    all_splice_structures = defaultdict(set)
    splice_struct_dict = defaultdict(list)
    for tx_info, exons in tx_dict.items():
        [chrom, strand, gene_id, tx_id] = tx_info
        exons.sort()
        splice_structure = []
        i = 0
        while (i + 1) < len(exons):
            introns = (exons[i][1] + 1, exons[i+1][0] - 1)
            splice_structure.append(introns)
            i = i + 1
        splice_structure = tuple(splice_structure)
        all_splice_structures[(chrom, strand)].add(splice_structure)
        if (chrom, strand) not in splice_struct_dict:
            splice_struct_dict[(chrom, strand)] = defaultdict(list)
        splice_struct_dict[(chrom, strand)][splice_structure].append(tx_info[2:])
    return all_splice_structures, splice_struct_dict

def process_tsv_splice_structures(tsv_file, min_int_count = 2):
    with gzip.open(tsv_file, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        all_splice_structures = defaultdict(set)
        for line in tbl:
            if not line[0].startswith('#'):
                [chrom, strand, introns] = line
                introns = introns.split(',')
                if len(introns) >= min_int_count:
                    new_introns = []
                    for intron in introns:
                        intron = re.sub('[\[\]()]', '', intron)
                        intron = tuple(map(int, intron.split('..')))
                        new_introns.append(intron)
                    all_splice_structures[(chrom, strand)].add(tuple(new_introns))
    return all_splice_structures

def write_filtered_gtf(gtf_in, gtf_out, all_feats, feats_dict):
    """
    Parse all_feats to get a list of feats to extract from gtf_in; using
    feats_dict as the dictionary. Extract all of those feats from gtf_in and
    write them to output gtf_out.
    """
    blessed_feats = defaultdict(set)
    for (chrom, strand), feats in all_feats.items():
        for feat in feats:
            tx_info = feats_dict[(chrom, strand)].get(feat, None)
            if tx_info:
                for tx in tx_info:
                    blessed_feats[(chrom, strand)].add(tx)
    with open(gtf_in, 'rt') as fi, open(gtf_out, 'wt') as fo:
        tblin = csv.reader(fi, delimiter = '\t')
        tblout = csv.writer(fo, delimiter = '\t', lineterminator = os.linesep, quotechar = '\xb6')
        for line in tblin:
            if len(line) != 9:
                tblout.writerow(line)
            else:
                [ chrom, source, feat_type,
                  start, stop, score,
                  strand, phase, attribs ] = line
                attribs = process_attribs(attribs)
                gene_id = attribs['gene_id']
                tx_id = attribs['transcript_id']
                if (gene_id, tx_id) in blessed_feats[(chrom, strand)]:
                    tblout.writerow(line)

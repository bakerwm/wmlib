#!/usr/bin/env python3

# convert gtf to bed6

# filter by feature type: 
# miRNA
# ncRNA
# nontranslated_CDS
# pre_miRNA
# protein_coding
# pseudogene
# rRNA
# snoRNA
# snRNA
# tRNA


import os
import sys
import re

def parse_desc(x, key='gene_id'):
    """
    Description in GTF:
    
    gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; transcript_name "4933401J01Rik-201"; transcript_source "havana"; transcript_biotype "TEC"; tag "basic"; transcript_support_level "NA";
    """
    if isinstance(x, str):
        d = {i.split()[0]:i.split()[1] for i in x.strip().split('; ')}
        out = d.get(key, None)
        if isinstance(out, str):
            out = out.replace('"', '')
    else:
        out = None
    return out


def gtf2bed(fh, feature='gene', gene_biotype='protein_coding'):
    """
    feature str
        gene|exon|CDS|transcript|five_prime_utr|three_prime_utr|start_codon|stop_codon|Selenocysteine
    gene_biotype str
        TEC|protein_coding|
    """
    for l in fh:
        p = l.strip().split('\t')
        # filtering
        if not len(p) == 9:
            continue
#         if 'mito' in p[0]:
#             continue
        if not p[2] == feature:
            continue
        # gene_name/gene_id
        name1 = parse_desc(p[8], 'gene_name')
        # name2 = parse_desc(p[8], 'gene_id')
        chr_name = p[0] if p[0].startswith('chr') else 'chr'+p[0]
        start = str(int(p[3]) - 1)
        bed = '\t'.join([chr_name, start, p[4], name1, '254', p[6]])
        yield bed


def main():
    if len(sys.argv) < 3:
        print('Usage: python gtf2bed <in.gtf> <out.bed>')
        sys.exit(1)
    gtf_in, bed_out = sys.argv[1:3]
    try:
        with open(gtf_in) as r, open(bed_out, 'wt') as w:
            for bed in gtf2bed(r):
                w.write(bed+'\n')
    except:
        print('failed read gtf: {}'.format(gtf_in))


if __name__ == '__main__':
    main()
# 

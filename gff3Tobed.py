#!/usr/bin/env python3
"""
Convert gff3 file to bed, gtf 
download from miRBase
"""

import os, sys, re


def gff3_reader(s):
    """
    read gff3 description
    ID=MIMAT0014959;Alias=MIMAT0014959;Name=mmu-miR-3113-5p;Derives_from=MI0014111
    return dict{}
    """
    p = s.strip().split('\t')
    if len(p) < 9:
        print(s)
    des = p[8].split(';')
    dd = {'chr': p[0],
        'start': p[3],
        'end': p[4],
        'strand': p[6],
        'feature': p[2],
        'source': p[1]}
    for d in des:
        key, val = d.split('=')
        dd[key] = val
    #
    return dd


def gff3_to_bed(d):
    """
    convert gff3 to bed
    input: dict
    output: string
    """
    chr = d['chr']
    start = str(int(d['start']) - 1)
    end = d['end']
    strand = d['strand']
    name = '%s,%s' % (d['Name'], d['ID'])
    bb = '\t'.join([chr, start, end, name, '100', strand])
    return bb


def gff3_to_gtf(d):
    """
    convert gff3 to gtf
    input: dict
    output: string
    """
    name = '%s,%s' % (d['Name'], d['ID'])
    des = 'gene_id "%s"; gene_name "%s"' % (name, name)
    g = [d['chr'], d['source'], d['feature'], d['start'],
        d['end'], '.', d['strand'], '.', des]
    return '\t'.join(g)


def bed_reader(s):
    """
    read bed line
    output: dict
    """
    p = s.strip().split('\t')
    if len(p) < 6:
        print(s)
    dd = {'chr': p[0],
        'start': p[1],
        'end': p[2],
        'name': p[3],
        'score': p[4],
        'strand': p[5]}
    return dd


def bed_to_gtf(d):
    """
    convert bed to gtf
    input: dict
    output: string
    """
    start = str(int(d['start']) + 1)
    des = 'gene_id "%s"; gene_name "%s"' % (d['name'], d['name'])
    g = [d['chr'], 'bed', 'exon', d['start'], d['end'], '.', 
        d['strand'], '.', des]
    return '\t'.join(g)


def gff3_converter(fn, fout, format='bed', feature='mature'):
    """
    read gff3 file, filt 'comment' lines, specific features
    name
    id
    """
    with open(fn, "rt") as fi, open(fout, 'wt') as fo:
        for line in fi:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if feature == 'mature':
                if not p[2] == 'miRNA':
                    continue
            else:
                if p[2] == 'miRNA':
                    continue
            d = gff3_reader(line)
            # convert to bed
            b = gff3_to_bed(d)
            # save to bed
            fo.write(b + '\n')
    print('save bed to file: %s' % fout)


def bed_converter(fn, fout, format='gtf'):
    """
    read bed file, 
    output: gtf
    """
    with open(fn, 'rt') as fi, open(fout, 'wt') as fo:
        for line in fi:
            d = bed_reader(line)
            g = bed_to_gtf(d)
            fo.write(g + '\n')
    print('save gtf to file: %s' % fout)


# convert to mature, precursor
if len(sys.argv) < 2 :
    sys.exit('usage: gff2bed.py in.gff3')

infile = sys.argv[1]
prefix, ext = os.path.splitext(os.path.basename(infile))

if ext == '.gff3':
    f_mature = prefix + '.mature.bed'
    f_precursor = prefix + '.precursor.bed'
    gff3_converter(infile, f_mature, feature='mature')
    gff3_converter(infile, f_precursor, feature='precursor')
elif ext == '.bed':
    f_gtf = prefix + '.gtf'
    f_precursor = prefix + '.gtf'
    bed_converter(infile, f_gtf)
else:
    sys.exit('unknown format: %s' % ext)

# convert to precursor    

#!/usr/bin/env python3

"""
Convert EMBL to GTF, GFF, Fasta ... formats

Date: 2020-07-09

Caution: 

Require modified version of biopython for EMBL file, 
Support: secondary_identifier

Add the following two lines after line-941 in `Scanner.py` file

if len(parts) == 3:
    consumer.dblink("%s:%s:%s" % (parts[0].strip(), parts[1].strip(), parts[2].strip())) # secondary_identifier

"""

import os
import sys
import re
from Bio import SeqIO # modified

# For transposon_elements, from FlyBase only

def get_dbxref(x, db=None):
    """
    Get the dbxrefs 
    dbxrefs format: {db}:{primary}:{secondary}
    """
    if db is None:
        ref = x.dbxrefs[0] if len(x.dbxrefs) > 0 else None
    else:
        ref = [i for i in x.dbxrefs if i.startswith(db)]
        ref = ref[0] if len(ref) > 0 else None

    # check ref
    if ref is None:
        parts = (None, None, None)
    else:
        parts = ref.split(':')
        if len(parts) < 3:
            parts.append(None)

    return parts


def embl2fa(input, output=None):
    """
    Convert EMBL to Fasta format
    dbxrefs: could be multiple
    """
    out = []
    for r in SeqIO.parse(input, 'embl'):
        xref = get_dbxref(r)[-1] # Dmel\\Idefix
        # only for Dmel
        if not xref.startswith('Dmel'):
            continue
        xref = xref.split('\\')[-1] # 
        out.append('>{}\n{}'.format(xref, str(r.seq)))

    if output is None:
        print('\n'.join(out))
    else:
        with open(output, 'wt') as w:
            write('\n'.join(out) + '\n')


def embl2tab(input, output=None):
    """
    Extract info
    id name flybase xref organism
    """
    out = []
    for r in SeqIO.parse(input, 'embl'):
        db, pri, sec = get_dbxref(r)
        org, xref = sec.split('\\') # organism, xref
        # out.append('\t'.join([r.id, r.name, db, pri, xref, org]))
        out.append('\t'.join([db, org, pri, xref, r.id, r.name]))

    if output is None:
        print('\n'.join(out))
    else:
        with open(output, 'wt') as w:
            write('\n'.join(out) + '\n')


def embl2gtf(input, output=None):
    """
    Extract info
    GTF:
    """
    out = []
    for r in SeqIO.parse(input, 'embl'):
        db, pri, sec = get_dbxref(r)
        org, xref = sec.split('\\') # organism, xref
        des = '; '.join([
            'gene_id \"{}\"'.format(xref),
            'gene_name \"{}\"'.format(xref),
            'transcript_id \"{}\"'.format(xref),
            'organism \"{}\"'.format(org),
            'dbxref \"{}\"'.format(pri)]) # !!!!
        gtf_gene = '\t'.join([
                xref,
                db,
                'gene',
                str(1),
                str(len(r.seq)),
                '.',
                '+',
                '.',
                des])
        gtf_exon = '\t'.join([
                xref,
                db,
                'exon',
                str(1),
                str(len(r.seq)),
                '.',
                '+',
                '.',
                des])
        out.append('{}\n{}'.format(gtf_gene, gtf_exon))

    if output is None:
        print('\n'.join(out))
    else:
        with open(output, 'wt') as w:
            write('\n'.join(out) + '\n')


def main():
    if len(sys.argv) < 3:
        sys.exit('read_embl.py [fa|gtf|text] input.embl > output')

    fmt = sys.argv[1].lower()
    embl = sys.argv[2]

    if fmt in ['fa', 'fasta']:
        embl2fa(embl)
    elif fmt in ['gtf']:
        embl2gtf(embl)
    elif fmt in ['tab', 'text']:
        embl2tab(embl)
    else:
        sys.exit('{} - unknown fmt, fa, gtf, tab expected')


if __name__ == '__main__':
    main()

#

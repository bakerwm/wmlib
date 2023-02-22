#!/usr/bin/env python3

import os
import sys
from xopen import xopen


if len(sys.argv) < 2:
	sys.exit('bed2gtf.py <in.bed>')

bed = sys.argv[1]

def bed2gtf(x, feature='gene'):
    """
    x : list
        bed record, BED6
    """
    if len(x) < 3:
        return None
    n1 = '{}:{}-{}'.format(x[0], x[1], x[2])
    name, strand = [x[3], x[5]] if len(x) > 5 else [n1, '+']
    s, e = x[1:3] # start, end
    s = int(s) + 1
    des = 'gene_id "{}"; gene_name "{}"'.format(name, name)
    gtf = [x[0], 'bed', feature, s, e, '.', strand, '.', des]
    return list(map(str, gtf))


with xopen(sys.argv[1]) as r:
    for l in r:
        g = bed2gtf(l.strip().split('\t'))
        if g:
            print('\t'.join(g))


# n = 0
# with open(bed) as fi:
#     for line in fi:
#         tab = line.strip().split('\t')
#         if len(tab) < 3:
#             continue
#         n += 1 #
#         # bed3 or bed6
#         if len(tab) > 3:
#             name = tab[3]
#             strand = tab[5]
#         else:
#             name = 'g{:06d}'.format(n)
#             strand = '+'
#         start = int(tab[1]) + 1
#         end = tab[2]
#         des = 'gene_id "{}"; gene_name "{}";'.format(name, name)
#         gtf = '\t'.join([
#             tab[0],
#             'bed',
#             'gene',
#             str(start),
#             end,
#             '.',
#             strand,
#             '.',
#             des])
#         print(gtf)

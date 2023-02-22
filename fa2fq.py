#!/usr/bin/env python3

"""
Convert fasta to fastq
quality I (40)
"""

import os
import sys
from xopen import xopen


def readfq(fh): # this is a generator
    """
    source: https://github.com/lh3/readfq/blob/master/readfq.py
    processing fastq file
    modified, return comment filed
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fh: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        [name, _, comment], seqs, last = last[1:].partition(" "), [], None
        for l in fh: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None, comment # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fh: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs), comment; # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None, comment # yield a fasta record instead
                break


def fa2fq(fa, fq):
    try:
        with xopen(fa) as r, xopen(fq, 'wt') as w:
            for s in readfq(r): # name,seq,qual,comment
                s2 = '@{}\n{}\n+\n{}'.format(s[0], s[1], 'I'*len(s[1]))
                w.write(s2+'\n')
    except:
        print('Could not read/write file: {}, {}'.format(fa, fq))


def main():
    if len(sys.argv) < 3:
        print('Usage: python fa2fq.py in.fa out.fq')
        sys.exit(1)
    fa2fq(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()

#
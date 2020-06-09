#!/usr/bin/env python3
"""
Tools for demx
1. check the top 10 barcode in each sample, saved in "comment" line
2. check barcode compatible for barcodes: (use a shiny-app)
"""

import os
import sys
import re
import fnmatch
from xopen import xopen
import pandas as pd


def listdir(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break

    return out


def listfile(path='.', pattern='*', full_name=True, recursive=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfile('./', '*.fq')
    """
    fn_list = listdir(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(f, pattern)]
    return fn_list


def readfq(fh): # this is a generator function
    """
    source: https://github.com/lh3/readfq/blob/master/readfq.py
    processing fastq file
    WM: Add comment filed to the end, 2020-04-12
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fh: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        tabs = last[1:].partition(" ")
        name, seqs, last = last[1:].partition(" ")[0], [], None
        comment = tabs[2] if len(tabs) > 2 else None # for comment field in name
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


def indexTable(x, topN=10):
    """Calculate the index table for fastq file"""
    ## extract the barcode
    d = {} # save barcode
    with xopen(x) as r:
        for name, seq, qual, comment in readfq(r):
            p = re.compile('[12]:[A-Z]:\d:(\w{8})')
            try:
                m = p.match(comment)
                bc = m.group(1)
            except:
                bc = comment.split(':')[-1]
            # out
            d[bc] = d.get(bc, 0) + 1

    ## convert to dataframe
    df = pd.DataFrame.from_dict(d, orient='index')
    df.columns = ['count']

    # sort by count
    df = df.sort_values(by='count', ascending=False)
    df = df.reset_index() # index to column

    # name
    df['id'] = re.sub('.gz$|.fq|.fastq', '', os.path.basename(x))

    # rename 
    df.columns = ['sequence', 'count', 'id']
    df = df.reset_index() # add index

    return df.head(topN)


if len(sys.argv) < 2:
    sys.exit('Usage: demxUtil.py <in.dir>')

indir = sys.argv[1]

fq_list = listfile(indir, "*fq.gz")
# fq_list = listdir(indir)

frames = [indexTable(f) for f in fq_list]
df = pd.concat(frames, axis=0, ignore_index=True)

print(df)

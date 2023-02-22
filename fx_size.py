#!/usr/bin/env python3

"""
read size distribution

Input; fastq, fasta, bam
"""

import os
import sys
import shutil
import pandas as pd
import numpy as np
import tempfile
import logging
import hiseq
from xopen import xopen

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def readfq(fh): # this is a generator function
    """
    source: https://github.com/lh3/readfq/blob/master/readfq.py
    processing fastq file
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


class FxSize(object):

    def __init__(self, fx, chunk_size=1000000, max_records=-1):
        self.fx = fx
        self.chunk_size = chunk_size
        self.max_records = max_records
        self.config()

        self.freq_table = self.fx_size()


    def config(self):
        # file exists
        if isinstance(self.fx, str):
            if not os.path.exists(self.fx):
                raise ValueError('fragsize.py, fx file not exists: {}'.format(self.fx))
        else:
            raise ValueError('fragsize.py, fx expect str, got {}'.format(type(self.fx).__name__))

        # file type
        self.fx_type = self.file_type(self.fx)
        self.is_fasta = self.fx_type == 'fasta'
        self.is_fastq = self.fx_type == 'fastq'
        if not self.fx_type in ['fasta', 'fastq']:
            raise ValueError('Not a fq|fa file: {}'.format(self.fx))

        # label
        fname = os.path.basename(self.fx)
        if fname.endswith('.gz'): # trim .gz
            fname = os.path.splitext(fname)[0]
        self.fname = os.path.splitext(fname)[0] # trim extension
        

    def file_type(self, fn, top_n=1000):
        """
        Check the file type by top 10000 rows:
        identify @ for fastq, > for fasta, * unknown
        """
        assert isinstance(fn, str)

        d = {}
        counter = 0
        with xopen(fn) as fh:
            for line in fh:
                counter += 1
                if counter > top_n:
                    break
                elif counter % 4 == 1: # for 1st line of fastq; and fa
                    base = line[0] # the first character
                    if base.lower() in 'acgtn': # sequence line in fasta
                        continue
                    d[base] = d.get(base, 0) + 1
                else:
                    continue

        ## percentage
        x = sorted(d.items(), key=lambda kv:kv[1], reverse=True)
        ## the top1 character
        x_top1 = x[0][0]
        x_top1_pct = x[0][1] / sum(d.values())

        ## check
        if x_top1 == '@':
            fx_type = 'fastq'
        elif x_top1 == '>':
            fx_type = 'fasta'
        else:
            fx_type = None

        ## if top1_pct < 90%
        if x_top1_pct < 0.9:
            fx_type = None

        return fx_type


    def cal_freq(self, x):
        """Calculate the frequency of list
        return dataframe

        index count
        """
        if isinstance(x, list):
            var, freq = np.unique(x, return_counts=True)
            df = pd.DataFrame(data=freq, index=var, columns=['count'])
        else:
            df = pd.DataFrame(columns=['count'])

        return df


    def fx_size(self):
        """
        size of fasta/q
        """
        flag = 0 # counter
        d = {}
        rs = [] # read size
        df = pd.DataFrame(columns=['count'])

        with xopen(self.fx, 'rt') as r:
            for _, seq, _, _ in readfq(r):
                rs.append(len(seq))
                flag += 1

                # maxmium
                if self.max_records > 0 and flag > self.max_records:
                    log.info('Stop at: {}'.format(flag))
                    break # stop

                # chunk
                if flag > 0 and flag % self.chunk_size == 0:
                    dfn = self.cal_freq(rs)
                    df = pd.concat([df, dfn], axis=1).sum(axis=1)
                    rs = [] # empty
                    log.info('{} : {} {}'.format('Processed', flag, self.fname))

        # last chunk
        if len(rs) > 0:
            dfn = self.cal_freq(rs)
            df = pd.concat([df, dfn], axis=1).sum(axis=1)
            rs = [] # empty
            log.info('{} : {} {}'.format('Processed', flag, self.fname))

        # convert to data.frame
        df = df.reset_index()
        df.columns = ['length', 'count']
        df['id'] = self.fname

        return df


    def distribution(self):
        """Basic statistics values
        value + freq

        mean, medium, mode, std, min, max, Q1, Q2, Q3
        """
        val = self.freq_table['length']
        freq = self.freq_table['count']
        inserts = np.repeat(val, freq)
        # inserts = np.repeat(self.freqTable['length'], self.freqTable['count'])

        # statistics
        q_mean = np.mean(inserts)
        q_median = np.median(inserts)
        q_median_dev = np.median(np.absolute(inserts - q_median))
        q_mode = val[np.argmax(freq)]
        q_std = np.std(inserts)
        q_min = np.min(inserts)
        q_max = np.max(inserts)
        q_qual = np.quantile(inserts, [0.25, 0.5, 0.75], axis=0)

        # core distribution
        s = np.array([q_mean, q_median, q_mode, q_std, q_min, q_max])
        s = np.append(s, q_qual)
        s = np.around(s, decimals=1)

        return s


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.csv',
            delete=False)
        return tmp.name


    def saveas(self, csv_file=None):
        """Save to file"""
        if csv_file is None:
            csv_file = self._tmp()

        log.info('saving to file: {}'.format(csv_file))

        try:
            # save table
            self.freq_table.to_csv(csv_file, index=False)

            # save statistics
            stat_file = csv_file + '.stat'
            da = self.distribution()
            pd.DataFrame(da).T.to_csv(stat_file, sep="\t", index=False, 
                header=['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1', 'Q2', 'Q3'])
        except:
            log.warning('failed saving file: {}'.format(csv_file))


    def plot(self, plot_file, x_min=0, x_max=1000, plot_type="bar"):
        """Generate freq table plot
        line plot, hist
        """
        # save to csv file
        csv_file = os.path.splitext(plot_file)[0] + '.csv'
        self.saveas(csv_file) # save csv, stat

        # save to plot
        # to-do: matplotlib function
        hiseq_dir = os.path.dirname(hiseq.__file__)
        fragPlotR = os.path.join(hiseq_dir, 'bin', 'qc.lendist.R')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            '{} {} {}'.format(fragPlotR, csv_file, plot_file),
            '{} {} {}'.format(plot_type, x_min, x_max)])
        print(cmd)

        try:
            os.system(cmd)
        except:
            log.warning('fragsize.py failed')


    def run(self):
        """Generate dataframe, pdf
        bam, str or list
        """
        pass


def main():
    if len(sys.argv) < 3:
        print('fx_stat.py <in.fq> <outdir> {<x_min:0> <x_max:1000>} ')
        sys.exit('args failed')

    fq, outdir = sys.argv[1:3]

    # output
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # cal fragsize
    f = FxSize(fq)

    # csv, plot
    plot_file = os.path.join(outdir, f.fname + '.lendist.pdf')
    f.plot(plot_file, x_min = 18, x_max = 40)


if __name__ == '__main__':
    main()












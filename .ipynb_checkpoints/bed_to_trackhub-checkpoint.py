#!/usr/bin/env python3
#-*- encoding: utf-8 -*- 


"""
Generate trackhub url for specific regions (BED)

mirror: http://genome.ucsc.edu
db: 
hubUrl:
pos:

example:
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1:1000-3000&hubUrl=<path_to_hub.txt> 

HubUrl(hub_url=hub_url, genome='dm6',
    position='chr2:10,000,000-12,000,000').trackhub_url()
"""


import os
import shutil
import argparse
import hiseq
from hiseq.get_trackhub.hub_url import HubUrl
from hiseq.utils.utils import log, url_to_link, run_shell_cmd
from hiseq.utils.file import (
    check_file, check_path, file_prefix, file_abspath, file_exists
)


def bed_to_trackhub(x, **kwargs):
    """
    load BED file, return the position: chr:star-end
    <generator>
    """
    # default args
    args = {
        'mirror': 'usa',
        'genome': None,
        'hub_url': None,
        'format': 'html',
        'append': False,
    }
    args.update(kwargs)
    # run
    try:
        with open(x) as r:
            for line in r:
                s = line.strip().split('\t')
                if len(s) >= 3:
                    p = '{chr}:{start}-{end}'.format(chr=s[0], start=s[1], end=s[2])
                    if len(s) > 3:
                        s[3] = os.path.basename(s[3]) # fix for narrowPeak
                        name = s[3]
                    else:
                        name = p
                    u = HubUrl(mirror=args['mirror'], 
                               genome=args['genome'], 
                               hub_url=args['hub_url'], 
                               position=p).trackhub_url()
                    k = url_to_link(u, name, args['format'])
                    if args['append']:
                        # k = ','.join([s[0], s[1], s[2], p, k])
                        k = ','.join(s + [p, k])
                    yield k                    
    except:
        log.error('failed reading file: {}'.format(x))

        
def report(**kwargs):
    """
    Generate report in html format
    
    <config>
    genome
    bed
    ...
    
    <table>
    chr,start,end,name,link
    """
    # output .csv, .html
    # default args
    args = {
        'input': None,
        'output': None,
        'mirror': 'usa',
        'genome': None,
        'hub_url': None,
        'format': 'html',
        'append': True,
    }
    args.update(kwargs)
    # check
    args['output'] = file_abspath(args['output']) #
    outdir = os.path.dirname(args['output'])
    if check_file(args['input']) and check_path(outdir):
        prefix = file_prefix(args['output'])
        out_html = os.path.join(outdir, prefix+'.table.html')
        # save to csv
        with open(args['output'], 'wt') as w:
            for p in bed_to_trackhub(args['input'], **args):
                w.write(p+'\n')
        # save to html
        pkg_dir = os.path.dirname(hiseq.__file__)
        report_R = os.path.join(pkg_dir, 'bin', 'csv_to_html.R')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            report_R,
            args['output'],
            out_html,
        ])
        # report_html
        if file_exists(out_html):
            log.info('report() skipped, file exists: {}'.format(out_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('report() failed')


def get_args():
    parser = argparse.ArgumentParser(
        prog='bed_to_trackhub',
        description='Generate url for specific position on trackhub',
        epilog='Example: \n\
               python trackhub_to_position.py')
    parser.add_argument('-i', '--input', required=True, default=None,
        help='Path to the BED file')
    parser.add_argument('-o', '--output', default=None,
        help='Path to the output file (.csv)')
    parser.add_argument('-m', '--mirror', default='usa', 
        help="The UCSC genome browser mirror, one of ['usa', 'euro', 'asia'] \
        or the URL of other mirror, default: [usa]")
    parser.add_argument('-g', '--genome', default=None, 
        help='The genome relase of the tracks, if [None], parsing this value \
        from the hub_url, in this way: hub_url -> hub.txt -> genomes.txt. \
        default: [None]')
    parser.add_argument('-l', '--hub-url', dest='hub_url', default=None,
        help='URL direct to the hub.txt file, that is open accessiable')
    parser.add_argument('-f', '--format', default='html', 
        choices=['html', 'markdown', 'url'],
        help='return the url in html/markdown/direct format, default: [html]')
    return parser


def main():
    args = vars(get_args().parse_args())
    report(**args)


if __name__ == '__main__':
    main()


#
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
import argparse
from hiseq.get_trackhub.hub_url import HubUrl
from hiseq.utils.utils import log


def bed_to_url(x, **kwargs):
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
    }
    args.update(kwargs)
    # run
    try:
        with open(x) as r:
            for line in r:
                s = line.strip().split('\t')
                if len(s) >= 3:
                    p = '{chr}:{start}-{end}'.format(chr=s[0], start=s[1], end=s[2])
                    name = s[3] if len(s) >= 4 else p # name filed
                    name = os.path.basename(name)  # fix for narrowPeak
                    u = HubUrl(mirror=args['mirror'], 
                               genome=args['genome'], 
                               hub_url=args['hub_url'], 
                               position=p).trackhub_url()
                    yield url_to_link(u, args['format'], name)
    except:
        log.error('unable to read file: {}'.format(x))

        
        
def url_to_link(x, format='markdown', n=None):
    """
    Convert url to link, in different format
    markdown, html
    
    markdown: [name](url) 
    html: <a href=url target='_blank'>name</a>
    """
    name = n if isinstance(n, str) else x
    if format == 'markdown':
        out = '[{name}]({url})'.format(name=name, url=x)
    elif format == 'html':
        out = "<a href={url} target='_blank'>{name}</a>".format(name=name, url=x)
    else:
        out = x
    return out
        

def get_args():
    parser = argparse.ArgumentParser(
        prog='trackhub_to_position',
        description='Generate url for specific position on trackhub',
        epilog='Example: \n\
               python trackhub_to_position.py')
    parser.add_argument('-i', '--input', required=True, default=None,
        help='Path to the BED file')
    parser.add_argument('-m', '--mirror', default='usa', 
        help="The UCSC genome browser mirror, one of ['usa', 'euro', 'asia'] \
        or the URL of other mirror, default: [usa]")
    parser.add_argument('-g', '--genome', default=None, 
        help='The genome relase of the tracks, if [None], parsing this value \
        from the hub_url, in this way: hub_url -> hub.txt -> genomes.txt. \
        default: [None]')
    parser.add_argument('-l', '--hub-url', dest='hub_url', default=None,
        help='URL direct to the hub.txt file, that is open accessiable')
    return parser
    
    
def main():
    args = vars(get_args().parse_args())
    for p in bed_to_url(args['input'], **args):
        print(p)


if __name__ == '__main__':
    main()


#
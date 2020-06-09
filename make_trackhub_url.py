#!/usr/bin/env python3

"""Share a track through URL

customTrack Params:

http://genome.ucsc.edu/cgi-bin/hgTracks?db=dm6&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrM%3A1%2D1000&hgsid=772022563_7rJ2OTnSjko4pS56BYAED3yJcEp7


base: http://genome.ucsc.edu/cgi-bin/

db={genome}

hubUrl={url-to-hub.txt}

position={chr:start-end}


# 1. search Trackhub directories
# 2. create hub.txt url
# 3. check url
# 4. generate share url


example-1: dir to hubUrl (hub.txt)

example-2: dir to hubShare (url)

example-3: specific position, dir to hubShare (url)

to-do:

support github, ...

"""

import os, sys, re
from pathlib import Path
import argparse
import subprocess
import logging


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)

log = logging.getLogger(__name__)


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub_url',
        description='Generate urls for sharing track_hubs',
        epilog='Example: ')
    parser.add_argument('--hubTxt', required=False, dest='hubTxt', default=None,
        help='The path to the track_hub hub.txt, if specified, ignore --source, --depth')
    parser.add_argument('--source', required=False, dest='source',
        help='The directory of track_hubs')
    parser.add_argument('--depth', required=False, metavar='depth', type=int,
        default=2, help='define the depth to search track_hub directories, default: [2]')
    parser.add_argument('--url-url', dest='url_url', 
        default='http://159.226.118.232/upload',
        help='The url direct to the track_hub directories; default: http://159.226.118.232/upload')
    parser.add_argument('--url-dir', dest='url_dir',
        default='/data/public/upload',
        help='The dir direct to the track_hub directories, default: [/data/public/upload]')
    # custome the track
    parser.add_argument('--hubCheck', required=False, action='store_true',
        help='if specified, check the hubUrl by command hubCheck, default: False')
    parser.add_argument('--position', required=False, default=None,
        help='The position to show in shareUrl. [chr:start-end]')
    parser.add_argument('--bed-in', required=False, metavar='bed_in',
        default=None, help='specify the position to show; format: chr:1-1000, require --hubTxt')
    parser.add_argument('--bed-out', required=False, metavar='bed_out',
        default=None, help='apend shareUrl to bed record, require --hubTxt')
    parser.add_argument('--log-level', default='INFO', dest='log_level',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')

    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def genome_reader(x):
    """Return the genome names from file *genomes.txt
    line:
    genome dm3
    ...
    """
    genomes = []
    with open(x) as fi:
        for line in fi:
            if line.startswith('genome'):
                tag, genome = line.strip().split()
                genomes.append(genome)
    return genomes


def list_dir(x, depth=2):
    """Search directories
    depth of the directory to search
    list all directories, with specific path depth
    """    
    dirs = []
    p = Path(x)
    # level = 1
    dirs = [p / x for x in p.iterdir() if x.is_dir()] # depth = 1
    depth -= 1 # 
    # sublevels
    if depth > 0:
        sub_dirs = []
        for i in dirs:
            sub_dirs.extend(list_dir(i, depth))
        # sub_dirs = [list_dir(i, depth) for i in dirs]
        dirs.extend(sub_dirs)

    return dirs


class HubDir(object):
    """Parse files in directory saving trackhub
    *genomes.txt
    *hub.txt
    {genome}/
    ...

    return:
    path-to-hub.txt
    genomes
    trackfiles
    """

    def __init__(self, path, url_url=None, url_dir=None, hubCheck=False,
        position=None):
        """required param
        path to the directory
        url_url: the website_path 
        url_dir: the directory to website_path
        hubCheck: use the command to check the hubUrl 
        """
        self.path = str(path)
        p = Path(str(path))
        self.files = [str(i) for i in p.iterdir() if i.is_file()]
        self.dirs = [str(i) for i in p.iterdir() if i.is_dir()]
        self.dirnames = [os.path.basename(i) for i in self.dirs]

        self.url_url = url_url
        self.url_dir = url_dir
        self.hubCheck = hubCheck
        self.position = position


    def hubTxt(self):
        """Search hub.txt file in path
        return absolute path
        """        
        hub_tmp = [i for i in self.files if re.search('hub.txt$', str(i))]
        if len(hub_tmp) == 0:
            # log.error('hub.txt not found in: {}'.format(self.path))
            hub_txt = None
        else:
            hub_txt = hub_tmp[0]

        return hub_txt


    def genomeTxt(self):
        """Search genomes.txt file in path
        return the absolute path
        """
        genome_tmp = [i for i in self.files if re.search('genomes.txt$', str(i))]
        if len(genome_tmp) == 0:
            # log.error('genomes.txt not found in: {}'.format(self.path))
            genome_txt = None
        else:
            genome_txt = genome_tmp[0]

        return genome_txt


    def genome(self):
        """Search the genomes.txt file
        record all the genomes 
        """
        genome_txt = self.genomeTxt()
        if genome_txt:
            genome = genome_reader(genome_txt)
            genome_check = [i for i in genome if i in self.dirnames]
            if not len(genome) == len(genome_check):
                genome = None
        else:
            genome = None

        return genome


    def is_hub(self):
        """Check if the directory is a hub
        *hub.txt
        *genomes.txt
        *genome
        """
        hub = self.hubTxt()
        genome_txt = self.genomeTxt()
        genome = self.genome()
        if hub and genome_txt and genome:
            tag = True
        else:
            tag = False

        return tag


    def is_hubUrl(self, x):
        """Validate the hubUrl by command : hubCheck
        check the file accessiable: hubUrl
        """
        hubCheck = which('hubCheck')

        if not hubCheck:
            log.error('command not found: hubCheck')
            return 0

        cmd = 'hubCheck -noTracks -checkSettings ' + x
        p1 = subprocess.run(cmd, shell=True, check=True, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        if p1.stderr.decode():
            log.error('hub-url error: {}'.format(x))
            return 0
        
        return x       


    def hubUrl(self):
        """Create url for track_hub on genome browser
        defautl: "track_hub"
        """
        hubTxt = self.hubTxt()

        # replace dir by url_url
        if self.url_url is None or self.url_dir is None or not self.url_dir in hubTxt:
            log.error('check url_url={}, url_dir={}'.format(self.url_url, selfurl_dir))
            url = None
        else:
            url = hubTxt.replace(self.url_dir, self.url_url)
            if self.hubCheck:
                url = self.is_hubUrl(url)

        return url


    def shareUrl(self, position=None,
        baseUrl='http://genome.ucsc.edu/cgi-bin/hgTracks'):
        """Create url for direct browser
        defail: 'web browser'
        """
        genomes = self.genome()
        genome = genomes[0]
        share = baseUrl + '?db={}'.format(genome)
        share += '&hubUrl={}'.format(self.hubUrl())
        position = self.position if position is None else position
        if not position is None:
            share += '&position={}'.format(position)

        return share


    def list_all(self):
        """List all fields
        save:
        genome:
        hubDir:
        hubUrl:
        shareUrl:
        """
        lines = []
        lines.append('{:>10}: {}'.format('genome', self.genome()[0]))
        lines.append('{:>10}: {}'.format('hubDir', self.path))
        lines.append('{:>10}: {}'.format('hubTxt', self.hubTxt()))
        lines.append('{:>10}: {}'.format('hubUrl', self.hubUrl()))
        lines.append('{:>10}: {}'.format('shareUrl', self.shareUrl()))

        return '\n'.join(lines)


def main():
    args = get_args()

    if args.hubTxt:
        hubDir = os.path.dirname(args.hubTxt)
        hub = HubDir(hubDir, url_url=args.url_url, url_dir=args.url_dir,
            hubCheck=args.hubCheck, position=args.position)
        
        # save share to BED
        if args.bed_in and args.bed_out:
            log.info('append shareUrl to BED file: {}'.format(args.bed_out))
            with open(args.bed_in) as f1, open(args.bed_out, 'wt') as f2:
                for line in f1:
                    line = line.strip()
                    chr, start, end = line.split('\t')[:3]
                    position = '{}:{}-{}'.format(chr, start, end)
                    shareUrl = hub.shareUrl()
                    f2.write(line + '\t' + shareUrl + '\n')

        hub_list = hub.list_all()
        print(hub_list + '\n')

    elif args.source:
        if args.bed_in or args.bed_out:
            log.warning('--bed-in, --bed-out, require --hubTxt')

        hubDirs = list_dir(args.source, args.depth)
        hubDirs = [i for i in hubDirs if HubDir(i).is_hub()]
        for x in hubDirs:
            hub = HubDir(x, url_url=args.url_url, url_dir=args.url_dir,
                hubCheck=args.hubCheck, position=args.position)
            hub_list = hub.list_all()
            print(hub_list + '\n')

    else:
        log.error('either --hubTxt or --source, required')



if __name__ == '__main__':
    main()





# def hubTxt2hubUrl(x, 
#     url_url='http://159.226.118.232/upload', 
#     url_dir='/data/public/upload'):
#     """Convert hubDir to hub.txt url
#     full-version
    
#     replace the url_dir by url_url in x
#     return new_url
#     """
#     # assert os.path.exists(x)

#     if url_dir in x:
#         url_new = x.replace(url_dir, url_url)
#     else:
#         url_new = None

#     return url_new

#     # # check URL
#     # if is_hubUrl(url_new):
#     #     return url_new
#     # else:
#     #     return 0


# def hubUrl2hubShare(x, genome, position=None):
#     """Share the hub, specific coordinates

#     ?db=genome
#     &hubUrl=url
#     &position=chr%3A1-100
#     """
#     baseurl = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
#     share = baseurl + '?db={}'.format(genome)
#     share += '&hubUrl={}'.format(x)
#     if position:
#         share += '&position={}'.format(position)

#     return share


# def is_hubDir(x):
#     """Check files in that directory
#     contains:
#     *genomes.txt
#     *hub.txt
#     genome/
#     """
#     p = Path(x)
    
#     if not os.path.exists(x):
#         return 0

#     # all files
#     x_files = [str(i) for i in p.iterdir() if i.is_file()]

#     # hub.txt file
#     hub_tmp = [i for i in x_files if re.search('hub.txt$', i)]

#     # read the genome.txt
#     genome_tmp = [i for i in x_files if re.search('genomes.txt$', i)]

#     # check-1
#     if len(hub_tmp) > 0 and len(genome_tmp) > 0:
#         hub_txt = str(p / hub_tmp[0]) # the first
#         genome_txt = str(p / genome_tmp[0]) # the first

#         # genome dir
#         genomes_dirs = []
#         if os.path.exists(genome_txt):
#             genomes = genome_reader(genome_txt)
#             genomes_dirs = [p / i for i in genomes]
#             genomes_dirs = [i for i in genomes_dirs if os.path.exists(p / i)]

#         # hub.txt, genomes.txt, {genome} dir
#         if len(genomes) == len(genomes_dirs):
#             return [genomes[0], hub_txt]

#     return 0


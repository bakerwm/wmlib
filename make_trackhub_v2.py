#!/usr/bin/env python
"""
This script is designed to create track_hub for bigWig + bigBed
composite tracks
make subgroups: strand (fwd/rev), replicate (if spcified), kind (signal, peak)


"""

import glob, os, sys, re, argparse
import subprocess, shlex, logging
import trackhub

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub',
        description='Generate trackhub for bigWig and bigBed files',
        epilog='Example: \
               python ~/work/wmlib/wmlib/make_trackhub.py --source clip_bigWig/ --hub PTB_CLIP \
               --genome hg19 --server yylab_tracks ')
    parser.add_argument('--source', required=True, dest='source',
        help='The directory of bigWig and bigBed files')
    parser.add_argument('--hub', required=True, metavar='hub_name',
        default=None, help='hub name (remove blanks)')
    parser.add_argument('--genome', required=True, metavar='GENOME', 
        default='hg19', help='genome for the trackhub')
    parser.add_argument('--short-label', default=None, dest='short_label',
        help='short label for the hub, default: [--hub]')
    parser.add_argument('--long-label', default=None, dest='long_label',
        help='long label for the hub, default: [--hub]')
    parser.add_argument('--user', default='Ming Wang',
        help='Who maintain the trackhub')
    parser.add_argument('--email', default='wangm08@hotmail.com',
        help='email for the hub')
    parser.add_argument('--server', default='yulab_tracks',
        help='Name of directory to save track files, default: [yulab_tracks]')
    parser.add_argument('--server-dir', dest='server_dir',
        default='/data/upload/ucsc_trackhub',
        help='Directory to save the track files, default: /data/upload/ucsc_trackhub')
    parser.add_argument('--server-url', default='http://159.226.118.232/open/',
        help='Server to deposite the data, default: http://159.226.118.232/open/')
    args=parser.parse_args()
    return args


##----------------------------------------------------------------------------##
## helper functions
# snippet is placed into public domain by
# anatoly techtonik <techtonik@gmail.com>
# http://stackoverflow.com/questions/8151300/ignore-case-in-glob-on-linux

import fnmatch
import os
import re

def findfiles(which, where='.'):
    """Returns list of filenames from `where` path matched by 'which'
    shell pattern. Matching is case-insensitive.
    # findfiles('*.ogg')
    """    
    # TODO: recursive param with walk() filtering
    rule = re.compile(fnmatch.translate(which), re.IGNORECASE)
    fn_names = [name for name in os.listdir(where) if rule.match(name)]
    return [os.path.join(where, f) for f in fn_names]


def listdir(x, recursive=False):
    """Return the files in path
    using os.walk()
    """
    file_list = []
    if recursive:
        for root, subdirs, files in os.walk(x):
            for d in subdirs:
                file_list.append(os.path.join(root, d))
            for f in files:
                file_list.append(os.path.join(root, f))
    else:
        file_list.extend(os.listdir(x))

    return file_list


def find_bw(x, bigbed=False):
    """Return the files in path"""
    assert isinstance(x, str)
    m1 = re.compile('\.bw$|\.bigwig', re.IGNORECASE)
    m2 = re.compile('\.bb$|\.bigbed', re.IGNORECASE)
    if bigbed:
        m=m2
    else:
        m=m1
    fn = [i for i in listdir(x, recursive=True) if m.search(os.path.basename(i))]
    return fn
    

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


def hub_url_checker(url):
    """Check the hubUrl by hubCheck"""
    hubCheck = which('hubCheck')
    tag = False
    if hubCheck is None:
        logging.info('hubCheck not found, please check the url manually.')
    else:
        p1 = subprocess.run(shlex.split('%s %s' % (hubCheck, url)),
            stderr=subprocess.PIPE)
        if p1.stderr.decode() == '':
            tag = True
        else:
            logging.error('hub-url error: %s' % url)
    return tag


def source_checker(path):
    """Check source directory, whether .bigWig or .bigBed files exists"""
    check_bw1 = glob.glob(os.path.join(path, '*.bigWig'))
    check_bw2 = glob.glob(os.path.join(path, '*.bw'))
    check_bg1 = glob.glob(os.path.join(path, "*.bigBed"))
    check_bg2 = glob.glob(os.path.join(path, "*.bb"))
    check_bw1.extend(check_bw2)

    return len(check_bw1) # require bigWig 


def pos_to_url(pos):
    """Convert position to url format
    Input: chr1:1-1000
    Output: &position=chr1%3A1-1000
    """
    assert isinstance(pos, str)
    pos = re.sub(',', '', p) # remove comma from string
    pos_url = None
    if re.compile('^chr\w+\:\d+\-\d+$').fullmatch(pos):
        chr, coords = pos.split(':')
        start, end = coords.split('-')
        pos_url = '&position=%s%3A%s-%s' % (chr, start, end)

    return pos_url


def hub_to_url(hub_txt, genome, position=None):
    """Create url for hub.txt
    hub_txt, the accessiable url of hub.txt file
    genome, UCSC supported name of genome, eg: hg19, dm3
    position, chr1:1-1000
    """
    assert isinstance(hub_txt, str)
    assert isinstance(genome, str)
    base_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
    track_url = '%s?db=%s&hubUrl=%s' % (base_url, genome, hub_txt)
    if position is None:
        return track_url
    else:
        pos_url = pos_to_url(position)
        return base_url + pos_url


def color_picker(fn):
    """Pick color based on the filename
    f:      forward, '#2E3440'
    r:      reverse, '#6DBFA7'
    other:  others, '#lalala'

    example: 
    rnaseq.fwd.bigWig
    
    to-do:
    assign multiple colors for tracks/files    
    """
    # Due to how code is extracted from the docs and run during tests, we
    # need to import again inside a function. You don't normally need this.
    import trackhub
    assert isinstance(fn, str)
    strand = get_strand(fn)

    colors = {
        'f': '#0000FE',
        'r': '#FE0000',
        'other': '#0000FE',
    }
    return trackhub.helpers.hex2rgb(colors[strand])


def build_subgroups():
    """Define the subgroups provide a way of tagging tracks."""
    subgroups = [
        # A subgroup for replicates
        trackhub.SubGroupDefinition(
            name='rep',
            label='rep',
            mapping={
                '0': 'merged',
                '1': 'rep1',
                '2': 'rep2',
                '3': 'rep3',
                '4': 'rep4',
                '5': 'rep5',
                '6': 'rep6',
            }
        ),

        # A subgroup for target
        trackhub.SubGroupDefinition(
            name='target',
            label='target',
            mapping={
                'piwi': 'shPIWI',
                'white': 'shWhite',
                'h3k9': 'H3K9me3',
                'other': 'Others',
            }
        ),

        # A subgroup for stage
        trackhub.SubGroupDefinition(
            name='stage',
            label='stage',
            mapping={
                '3h': '3h',
                '6h': '6h',
                '9h': '9h',
                '12h': '12h',
                'l1': 'L1',
                'l2': 'L2',
                'l3': 'L3',
                'other': 'Others',
            }
        ),

        # Forward/Reverse strand
        trackhub.SubGroupDefinition(
            name='strand',
            label='strand',
            mapping={
                'f': 'Forward',
                'r': 'Reverse',
            }
        ),

        # Turning on/off the signal or regions in bulk.
        trackhub.SubGroupDefinition(
            name='kind',
            label='kind',
            mapping={
                'signal': 'signal',
                'peak': 'peak',
            }
        ),
    ]
    return subgroups


def get_strand(fn):
    """Determine the strand of file based on filename"""
    assert isinstance(fn, str)
    fn_name = os.path.splitext(fn)[0].lower()
    r1 = re.compile('fwd$|watson$|plus$', re.IGNORECASE)
    r2 = re.compile('rev$|crick$|minus$', re.IGNORECASE)

    if r1.search(fn_name):
        strand = 'f'
    elif r2.search(fn_name):
        strand = 'r'
    else:
        strand = 'other'

    return strand


def get_target(fn):
    """Get the target of the file
    shPIWI or shWhite
    """
    assert isinstance(fn, str)
    fn_name = os.path.basename(fn).lower()
    if fn_name.startswith('h3k9'):
        target = 'h3k9'
    elif fn_name.startswith('shpiwi'):
        target = 'piwi'
    elif fn_name.startswith('shwhite'):
        target = 'white'
    else:
        target = 'other'

    return target


def get_stage(fn):
    """Get the stage of the file
    3h, 6h, 12h, L1, L2, L3
    """
    assert isinstance(fn, str)
    fn_name = os.path.basename(fn).lower()
    m = re.compile('\d+h_rep|\w\d_rep').search(fn_name)
    if m:
        stage = fn_name[m.start():m.end()]
        stage = re.sub(r'_rep', '', stage)
    else:
        stage = 'other'

    return stage


def get_rep(fn):
    """Get the replicate name of the file
    rep1, rep2, ...
    """
    assert isinstance(fn, str)
    fn_name = os.path.splitext(fn)[0].lower()
    m = re.compile('_rep[0-9]').search(fn_name)
    if m:
        rep = fn_name[m.start():m.end()]
        rep = re.sub('_rep', '', rep) # the number
    else:
        rep = 0

    return rep


def get_subgroup(fn):
    """Determine the subgroups by filename
    This functions figures out subgroups based on the number in the
    filename.  Subgroups provided to the Track() constructor is
    a dictionary where keys are `rep` attributes from the subgroups added
    to the composite above, and values are keys of the `mapping` attribute
    of that same subgroup.

    Might be easier to cross-reference with the subgroups above, but an
    example return value from this function would be:
    """
    assert isinstance(fn, str)
    # 1. strand
    strand = get_strand(fn)

    # 2. target
    target = get_target(fn)

    # 3. stage
    stage = get_stage(fn)

    # 4. kind, peak/signal
    fn_ext = os.path.splitext(fn)[1].lower()
    if fn_ext == '.bw' or fn_ext == '.bigwig':
        kind = 'signal'
    elif fn_ext == '.bb' or fn_ext == '.bigbed':
        kind = 'peak'
    else:
        kind = 'peak'
    
    # 5. replicate
    rep = get_rep(fn)

    ## construct subgroup
    track_subgroup = {
        'strand': strand,
        'target': target,
        'stage': stage,
        'kind': kind,
        'rep': rep,
    }

    return track_subgroup


def create_composite(trackdb, short_label, bigwig_files, 
    bigbed_files=None, bedgraph_files=None):
    """ Create composite tracks, bigWig, bigBed,
    Add CompositeTrack to trackdb
    """
    assert isinstance(trackdb, trackhub.trackdb.TrackDb)
    assert isinstance(short_label, str)
    assert isinstance(bigwig_files, list)

    # Create the composite track
    composite = trackhub.CompositeTrack(
        name='composite',
        short_label=short_label,
        long_label=short_label,
        dimensions='dimX=stage dimY=target dimA=kind dimB=rep dimC=strand',
        filterComposite='dimA',
        sortOrder='stage=+ target=+ kind=-',
        tracktype='bigWig',
        visibility='full'
    )

    # Add subgroups to the composite track
    composite.add_subgroups(build_subgroups())

    # Add the composite track to the trackDb
    trackdb.add_tracks(composite)

    # CompositeTracks compose different ViewTracks.
    signal_view = trackhub.ViewTrack(
        name='signalviewtrack',
        short_label='Signal',
        view='signal',
        visibility='full',
        tracktype='bigWig')

    # signal_view2 = trackhub.ViewTrack(
    #     name='singalviewtrack2',
    #     short_label='Signal2',
    #     view='signal',
    #     visibility='full',
    #     tracktype='bedGraph')

    regions_view = trackhub.ViewTrack(
        name='regionsviewtrack',
        short_label='Regions',
        view='regions',
        visibility='dense',
        tracktype='bigBed')

    # for bigWig files
    # signal view
    composite.add_view(signal_view)
    for bigwig in bigwig_files:
        track=trackhub.Track(
            name=trackhub.helpers.sanitize(os.path.basename(bigwig)),
            source=bigwig,
            visibility='full',
            tracktype='bigWig',
            viewLimits='0:15',
            maxHeightPixels='8:40:128',
            subgroups=get_subgroup(bigwig),
            color=color_picker(bigwig))

        # Note that we add the track to the *view* rather than the trackDb as
        # we did in the README example.
        signal_view.add_tracks(track)


    # # for bedgraph
    # if bedgraph_files is None:
    #     pass
    # elif len(bedgraph_files) > 0:
    #     composite.add_view(signal_view2)
    #     for bg in bedgraph_files:
    #         track = trackhub.Track(
    #             name=trackhub.helpers.sanitize(os.path.basename(bg)),
    #             source=bg,
    #             visibility='full',
    #             subgroups=get_subgroup(bg),
    #             color=color_picker(bg),
    #             tracktype='bedGraph',
    #             viewLimits='-2:2',
    #             maxHeightPixels='8:40:128')
    #         signal_view2.add_tracks(track)
    # else:
    #     pass


    # for bigBed files
    # region view
    # bigbed_files = glob.glob(os.path.join(data_path, '*.bigBed'))
    if bigbed_files is None:
        pass
    elif len(bigbed_files) > 0 :
        composite.add_view(regions_view)
        for bigbed in bigbed_files:
            track = trackhub.Track(
                name=trackhub.helpers.sanitize(os.path.basename(bigbed)),
                source=bigbed,
                visibility='full',
                subgroups=get_subgroup(bigbed),
                color=color_picker(bigbed),
                tracktype='bigBed')
            regions_view.add_tracks(track)
    else:
        pass

    return composite # CompositeTrack


def get_bigwig(x):
    """List bigwig files from x path"""
    


def main():
    args = get_args()

    # remove spaces for hub_name
    hub_name = re.sub('\s+', '', args.hub)

    # check labels for hub
    if args.short_label is None:
        args.short_label = hub_name
    if args.long_label is None:
        args.long_label = args.long_label

    # check server name
    if args.server is None:
        args.server = hub_name
    
    # server dir, save files
    server_path = os.path.join(args.server_dir, args.server, hub_name)
    if not os.path.exists(server_path):
        try:
            os.makedirs(server_path, 0o711)
        except IOError:
            print('cannot create directory: ' + server_path)

    # get files
    bw = find_bw(args.source)
    bb = find_bw(args.source, bigbed=True)

    #################
    # Create tracks #
    #################
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=args.short_label,
        long_label=args.long_label,
        genome=args.genome,
        email=args.email)

    #########################
    # Create CompositeTrack #
    #########################
    composite = create_composite(
        trackdb=trackdb, 
        short_label=args.short_label, 
        bigwig_files=bw, 
        bigbed_files=bb)

    trackhub.upload.upload_hub(hub=hub, host='localhost', 
        remote_dir=server_path)

    # hub_txt
    # server:     hub_name
    # server_dir: /data/ucsc_gb/trackhub_open/
    # server_url: http://159.226.118.232/open/
    # 
    # server_path: /data/ucsc_gb/trackhub_open/<server>/
    # server_url:  http://159.226.118.232/open/<server>/
    # hub_txt_url: http://159.226.118.232/open/<server>/<hub_name>.hub.txt
    hub_txt_url = os.path.join(args.server_url, args.server, hub_name, hub_name + '.hub.txt')
    hub_url = hub_to_url(hub_txt_url, args.genome)


    if hub_url_checker(hub_txt_url) is True:
        status = 'ok'
    else:
        status = 'failed'
    print('%10s : %s' % (status, hub_txt_url))
    print('%10s : %s' % (status, hub_url))


    # Example uploading to a web server (not run):
    if 0:
        trackhub.upload.upload_hub(
            hub=hub, host='example.com', user='username',
            remote_dir='/var/www/example_hub')

if __name__ == '__main__':
    main()


## EOF

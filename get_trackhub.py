#!/usr/bin/env python
"""
This script is designed to create track_hub for bigWig + bigBed
composite tracks
make subgroups: strand (fwd/rev), replicate (if spcified), kind (signal, peak)
"""

import glob, os, re, argparse
import trackhub

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'get_trackhub',
        description = 'Generate trackhub files',
        epilog = 'Example: python get_trackhub.py --data_dir data/ --hub demo --genome hg19 --short_label ptb_clip')
    parser.add_argument('--data_dir', required = True,
        help = 'The directory of bigWig and bigBed files')
    parser.add_argument('--hub', required = True, metavar = 'hub_name',
        default = None, help = 'hub name (remove blanks)')
    parser.add_argument('--genome', required = True, metavar = 'GENOME', 
        default = 'hg19', help = 'genome for the trackhub')
    parser.add_argument('--short_label', required = False,
        default = None, help = 'short label for the hub')
    parser.add_argument('--long_label', required = False,
        default = None, help = 'long label for the hub')
    parser.add_argument('--user', default = 'Ming Wang',
        help = 'Who maintain the trackhub')
    parser.add_argument('--email', default = 'wangm08@hotmail.com',
        help = 'email for the hub')
    parser.add_argument('--server', default = '/home/data/ucsc_gb/track_hub/',
        help = 'Server to save the files, default: /home/data/ucsc_gb/track_hub/')
    args = parser.parse_args()
    return args


##----------------------------------------------------------------------------##
## helper functions
def hub_to_url(hub_txt, genome, position=None):
    """
    create url for hub.txt
    position: chr1:1-1000 convert to url encode
    'chr1%3A1-1000'
    """
    base_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
    db = '?db=' + genome
    hubUrl = '&hubUrl=' + hub_txt
    if position:
        pos = '&position=' + re.sub(r':', r'%3A', position)
        return(base_url + db + hubUrl + pos)
    else:
        return(base_url + db + hubUrl)


def get_subgroups():
    # Subgroups provide a way of tagging tracks.
    subgroups = [

        # A subgroup to select the replicate number
        trackhub.SubGroupDefinition(
            name = 'rep',
            label = 'replicate',
            mapping = {
                '0': 'merged',
                '1': 'rep1',
                '2': 'rep2',
                '3': 'rep3',
                '4': 'rep4',
            }
        ),

        # Forward/Reverse strand
        trackhub.SubGroupDefinition(
            name = 'strand',
            label = 'strand',
            mapping = {
                'f': 'Forward',
                'r': 'Reverse',
            }
        ),

        # Turning on/off the signal or regions in bulk.
        trackhub.SubGroupDefinition(
            name = 'kind',
            label = 'kind',
            mapping = {
                'signal': 'signal',
                'peak': 'peak',
            }
        ),
    ]
    return subgroups


def subgroups_from_filename(fn):
    """
    This functions figures out subgroups based on the number in the
    filename.  Subgroups provided to the Track() constructor is
    a dictionary where keys are `rep` attributes from the subgroups added
    to the composite above, and values are keys of the `mapping` attribute
    of that same subgroup.

    Might be easier to cross-reference with the subgroups above, but an
    example return value from this function would be:
    """
    n = os.path.basename(fn).split('.')[0].split('_')
    strand = 'fwd' if len(n) < 2 else n[-2]
    # strand = os.path.basename(fn).split('.')[0].split('_')[-2]
    if strand == 'rev' or 'rev' in os.path.basename(fn):
        strand = 'r'
    else:
        strand = 'f'
    if fn.endswith('bw') or fn.endswith('bigWig'):
        kind = 'signal'
    else:
        kind = 'peak'
    # rep = os.path.basename(fn).split('.')[0].split('_')[-3]
    rep = '0' if len(n) < 3 else n[-3]
    if rep.startswith('rep'):
        rep = re.sub('rep', '', rep) # replicate
    else:
        rep = 0 # merged
    track_subgroup = {
        'rep': rep,
        'strand': strand,
        'kind': kind,
    }
    return track_subgroup


def color_from_filename(fn):
    """
    Figure out a nice color for a track, depending on its filename.
    """
    # Due to how code is extracted from the docs and run during tests, we
    # need to import again inside a function. You don't normally need this.
    import trackhub
    n = os.path.basename(fn).split('.')[0].split('_')
    strand = 'fwd' if len(n) < 2 else n[-2]
    # strand = os.path.basename(fn).split('.')[0].split('_')[-2]
    strand = strand if strand in ['fwd', 'rev'] else 'other'
    colors = {
        'fwd': '#2E3440',
        'rev': '#6DBFA7',
        'other': '#1a1a1a',
    }
    return trackhub.helpers.hex2rgb(colors[strand])


def create_composite(trackdb, track_name, comp_label, data_path):
    ##------------------------------------------------------------------------##
    ## compositeTrack
    # Create the composite track
    composite = trackhub.CompositeTrack(
        name = track_name,
        short_label = comp_label,

        # The available options for dimensions are the `name` attributes of
        # each subgroup. Start with dimX and dimY (which become axes of the
        # checkbox matrix to select tracks), and then dimA, dimB, etc.
        dimensions = 'dimX=strand dimY=rep dimA=kind',

        # Enable a drop-down box under the checkbox matrix
        filterComposite = 'dimA',

        # `name` attributes of each subgroup
        sortOrder = 'rep=+ kind=-',
        tracktype = 'bigWig',
        visibility = 'full'
    )

    # Add subgroups to the composite track
    composite.add_subgroups(get_subgroups())

    # Add the composite track to the trackDb
    trackdb.add_tracks(composite)

    # CompositeTracks compose different ViewTracks.
    # one for signal in bigWig, another one for bigBed regions.
    signal_view = trackhub.ViewTrack(
        name = track_name + 'signal',
        short_label = 'Signal',
        view = 'signal',
        visibility = 'full',
        tracktype = 'bigWig'
        )

    regions_view = trackhub.ViewTrack(
        name = track_name + 'view',
        short_label = 'Regions',
        view = 'regions',
        visibility = 'dense',
        tracktype = 'bigWig')

    # for bigWig files
    composite.add_view(signal_view)
    for bigwig in glob.glob(os.path.join(data_path, '*.bigWig')):
        track = trackhub.Track(
            name=trackhub.helpers.sanitize(os.path.basename(bigwig)),
            source=bigwig,
            visibility='full',
            tracktype='bigWig',
            viewLimits='-2:2',
            maxHeightPixels='8:40:128',
            subgroups=subgroups_from_filename(bigwig),
            color=color_from_filename(bigwig)
        )

        # Note that we add the track to the *view* rather than the trackDb as
        # we did in the README example.
        signal_view.add_tracks(track)

    # for bigBed files
    bigbed_files = glob.glob(os.path.join(data_path, '*.bigBed'))
    if len(bigbed_files) > 0 :
        composite.add_view(regions_view)
        for bigbed in bigbed_files:
            track = trackhub.Track(
                name=trackhub.helpers.sanitize(os.path.basename(bigbed)),
                source=bigbed,
                visibility='full',
                subgroups=subgroups_from_filename(bigbed),
                color=color_from_filename(bigbed),
                tracktype='bigBed'
            )
            regions_view.add_tracks(track)



def dir_checker(path):
    """
    check bigWig and bigBed files in path
    """
    check_bw = glob.glob(os.path.join(path, '*.bigWig'))
    check_bg = glob.glob(os.path.join(path, "*.bigBed"))
    # if len(check_bw) * len(check_bg) == 0:
    if len(check_bw) == 0:
        return False
    else:
        return True

def main():
    args = get_args()

    # remove spaces for hub_name
    hub_name = re.sub('\s+', '', args.hub)
    genome = args.genome
    hub_user = args.user
    hub_email = args.email

    # check labels for hub
    hub_short_label = args.short_label if args.short_label else args.hub
    hub_long_label = args.long_label if args.long_label else hub_short_label
    
    # data_dir
    data_dir = args.data_dir

    # upload dir
    upload_dir = os.path.join(args.server, hub_name)
    try:
        os.makedirs(upload_dir, 0o711)
    except IOError:
        print('cannot create directory: ' + upload_dir)

    # subfolders for sampels
    # each subfolders should contain *.bigWig and *.bigBed files
    if dir_checker(data_dir):
        sub_dirs = [data_dir]
    else:
        # sub_dirs in data_dir
        sub_dirs = [d for d in glob.glob(os.path.join(data_dir, '*')) 
                    if os.path.isdir(d)]
        sub_dirs = [d for d in sub_dirs if dir_checker(d)]
    assert len(sub_dirs) > 0

    # Initialize the components of a track hub
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name = hub_name,
        short_label = hub_short_label,
        long_label = hub_long_label,
        genome = genome,
        email = hub_email)

    # create composite tracks
    for d in sub_dirs:
        create_composite(trackdb, os.path.basename(d) + '_composite', 
                         os.path.basename(d), d)

    # upload data to remote_dir
    print(upload_dir)
    trackhub.upload.upload_hub(hub = hub, host = 'localhost', 
        remote_dir = upload_dir)

    # Example uploading to a web server (not run):
    if 0:
        trackhub.upload.upload_hub(
            hub=hub, host='example.com', user='username',
            remote_dir='/var/www/example_hub')

if __name__ == '__main__':
    main()


## EOF

#!/usr/bin/env python3 

"""Create tracks using pyGenomeTracks

Template for tracks.ini file

[test bigwig lines]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 2
type = line
title = orientation = inverted; show_data_range = false
orientation = inverted
show_data_range = false
max_value = 50

[spacer]

[test bigwig min]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
type = line
summary_method = min
max_value = 150
min_value = -25
overlay_previous = share-y
number_of_bins = 300

[spacer]

[x-axis]


# command
$ pyGenomeTracks  --tracks tracks.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o master_bigwig.png



##################
## targets   
1. bigWig + peak (BED/narrowPeak) + (genes: GTF)
2. bigWig + (genes: GTF)  
3. bigWIg (overlay, ...)
"""

import os
import sys
import pathlib
import numbers
import pysam
import pyBigWig
from hiseq.utils.helper import *



################################################################################
## functions for pipe() project
################################################################################
def get_n_map(x):
    """The number of mappped reads for project
    
    Parameters
    ----------
    x : str
        Path to the pipe() project directory
        
    group : str
        The group of reads, choose from ['map', 'total', 'clean', 'collapse',
        'smRNA', 'miRNA', 'te', 'piRC', 'genome', 'unmap'], default: `map`
    
    11.stat/fx_stat.reads.csv
    """
    # stat, reads
    r_stat = os.path.join(x, '11.stat', 'fx_stat.reads.csv')  
    try:
        df = pd.read_csv(r_stat)
        n = df['all'][1] - df['all'][11]  # 1.clean - 11.unmap
    except:
        print('failed to read file: {}'.format(r_stat))
        n = 0
    return n


def bam_to_bw(bam, bw, scale=1, gsize=0, threads=8):
    """Convert bam to bigWig, using bamCoverage (deeptools)
    Parameters
    ----------
    bam : str
        path to the bam file
        
    bw : str
        path to the bigWig file
        
    scale : float
        the scale factor, [0-1], default: 1
    """
    # check bam exists
    # return the gsize
    ref_size = gsize if gsize > 0 else sum(pysam.AlignmentFile(bam).header.lengths)
    # cmd
    cmd_sh = os.path.join(os.path.dirname(bam), os.path.basename(bam)+'.bam2bw.sh')
    cmd_log = os.path.join(os.path.dirname(bam), os.path.basename(bam)+'.bam2bw.log')
    cmd = ' '.join([
        'bamCoverage', 
        '-b {}'.format(bam),
        '-o {}'.format(bw),
        '--scaleFactor {}'.format(scale),
        '--effectiveGenomeSize {}'.format(ref_size),
        '-p {}'.format(threads),
        '--binSize 1 --normalizeUsing None',
        '2> {}'.format(cmd_log),
    ])
    with open(cmd_sh, 'wt') as w:
        w.write(cmd+'\n')
    # run
    if os.path.exists(bw):
        log.warning('get_te_bw() skipped, file exists: {}'.format(bw))
    else:
        os.system(cmd)
    return bw
    

def get_te_bam(x, index='te', unique='unique'):
    """Extract the TE, unique bam
    
    Parameters
    ----------
    x : str
        path to the pipe() project
    
    index : str
        the index name, ['te', 'piRC', 'genome'], default: 'te'
        
    unique : str
        unique mappping reads or not, ['unique', 'multi', 'both'], default: 'unique'
    """
    if not isinstance(x, str):
        sys.exit('x expect str, got {}'.foramt(type(x).__name__))
    config = os.path.join(x, 'config', 'config.toml')
    if not file_exists(config):
        sys.exit('x expect pip() dir, config not found: {}'.format(config))
    p = Config().load(config)
    # link to files
    index_dict = {
        'te': '07.te',
        'piRC': '08.piRNA_cluster',
        'piRNA_cluster': '08.piRNA_cluster',
        'genome': '09.genome',
    }
    index_name = index_dict.get(index, '07.te')    
    return os.path.join(x, index_name, unique, p['smp_name']+'.bam') # name

    
    
def get_te_bw(x, strandness=False):
    """
    Extract strandness bigWig files: (.fwd, .rev)
    
    Parameters
    ----------
    x : str
        path to the project_dir (pipe)

    # mapped reads: clean_data - unmap
    x/07.te/unique/*.bigwig
    """
    te_bam = get_te_bam(x, 'te', 'unique')
    prefix = os.path.splitext(te_bam)[0]
    # scale
    n_map = get_n_map(x)
    scale = 1e6/n_map if n_map > 0 else 1
    # TE
    te_bw = prefix+'.bigwig'
    if strandness:
        # forward
        te_bam_fwd = prefix+'.fwd.bam'
        te_bw_fwd = prefix+'.fwd.bigwig'
        if not file_exists(te_bam_fwd):
            pysam.view('-bhS', '-F', '16', '-o', te_bam_fwd, te_bam, catch_stdout=False)
            pysam.index(te_bam_fwd)
        bam_to_bw(te_bam_fwd, te_bw_fwd, scale=scale)
        # reverse
        te_bam_rev = prefix+'.rev.bam'
        te_bw_rev = prefix+'.rev.bigwig'
        if not file_exists(te_bam_rev):
            pysam.view('-bhS', '-f', '16', '-o', te_bam_rev, te_bam, catch_stdout=False)
            pysam.index(te_bam_rev)
        bam_to_bw(te_bam_rev, te_bw_rev, scale=scale)
        out = [te_bw_fwd, te_bw_rev]
    else:
        bam_to_bw(te_bam, te_bw, scale=scale)
        out = [te_bw, None]
    return out



################################################################################
## functions for pyGenomeTracks
################################################################################
def load_regions_from_bam(x):
    """Extract the regions for whole chromosome
    x is the bam file
    chr1:1-1000    
    """
    try:
        h = pysam.AlignmentFile(x).header
        out = {i:(1, k) for i,k in zip(h.references, h.lengths)}
    except:
        print('load_regions_from_bam() failed, {}'.format(x))
        out = None
    return out




def load_regions_from_bw(x):
    """Extract the chr 
    
    Parameters
    ----------
    x : str
        The bigWig file
    
    chr: str
        The name of the chromosome, if `None`, return the whole genome
    """
    try:
        bw = pyBigWig.open(x)
        d = bw.chroms()
        out = {k:(1, v) for k,v in d.items()}
    except:
        print('load_regions_from_bw() failed, {}'.format(x))
        out = None
    return out



def load_regions_from_bed(x):
    """Load regions from BED file
    BED3, BED6, ...
    
    Parameters
    ----------
    x : str
        BED file
    """
    try:
        d = {}
        for bed in pybedtools.BedTool(x):
            d[bed.chrom] = (int(bed.start)+1, int(bed.end))
    except:
        print('load_regions_from_bed() failed, {}'.foramt(x))
        d = None
    return d
            
    
    
    
    
    
def bw_to_config_f1(x, **kwargs):
    """Generate the config.ini file for files: bigWig

    Parameters
    ----------
    x : str
        path to the bigwig file

    Keyword arguments
    -----------------
        color=None, 
        min_value=0, 
        max_value='auto', 
        nbins=700

    Example output
    --------------
    [track-01]
    file = 
    title = 
    height = 2
    color = 
    min_value = auto
    max_value = auto
    number_of_bins = 700
    file_type = bigwig
    """
    # arguments
    xcolor = kwargs.get('color', '#333333') # default: grey
    min_value = kwargs.get('min_value', 0)
    max_value = kwargs.get('max_value', 'auto')
    nbins = kwargs.get('nbins', 700)
    if os.path.exists(x):
        xname, xtype = os.path.splitext(os.path.basename(x))
        xtype = xtype.lstrip('.')
        return '\n'.join([
            '[track-{}]'.format(xname),
            'file = {}'.format(x),
            'title = {}'.format(xname),
            'height = {}'.format(2),
            'color = {}'.format(xcolor),
            'min_value = {}'.format(min_value),
            'max_value = {}'.format(max_value),
            'number_of_bins = {}'.format(nbins),
            'file_type = {}'.format(xtype.lower()),
        ])
        

        
def bw_to_config(x, chr, **kwargs):
    """Generate the config.ini and run.sh for the region
    
    Parameters
    ----------
    x : list
        list of bigWig files

    chr : str
        specify the chr name

    Keyword arguments:
        outdir: Directory to save the config.ini and plot
                if `None`, set the `cwd()`.
        start:  Starting position
        end:    Ending position
        y_min:  The minimum value on y-axis, default: 0
        y_max:  The maximum value on y-axis, default: 'auto'
        colors: Set colors for the tracks/bigWig files, 
                if `auto`, use colors by `fishualize` package
                from R, or set a list of colors by name, 
                `['red', 'green', 'blue']`, the number should
                match the bigWig files
        config: The name of `ini` file, defualt: `outdir/config.ini`
        plot_file: The name of the plot file, 'chr_start_end.png'
        
    Example
    [track-01]
    ...
    [spacer]
    [x-axis]
    """
    if isinstance(x, list):
        x = [b for b in x if file_exists(b)] # remove None
    else:
        print('x expect list of bigWig files, got {}'.format(
            type(x).__name__))
        return None
    if len(x) == 0:
        print('x, list of bigWig files, not found')
        return None
    x = [file_abspath(i) for i in x]
    ###########################################################################
    ## arguments
    if not isinstance(chr, str):
        print('chr expect str, got {}'.format(type(chr).__name__))
        return None
    ## outdir
    outdir = kwargs.get('outdir', None)
    if not isinstance(outdir, str):
        outdir = str(pathlib.Path.cwd())
    outdir = file_abspath(outdir)
    project_dir = os.path.join(outdir, chr)
    check_path(project_dir)
    ## start/end
    start = kwargs.get('start', 1)
    end = kwargs.get('end', 1000)
    ## colors
    i_colors = kwargs.get('colors', None)
    if isinstance(i_colors, str): 
        x_colors = i_colors * len(x)
    elif isinstance(i_colors, list):
        if len(i_colors) >= len(x):
            x_colors = i_colors[:len(x)]
        else:
            x_colors = fish_colors(len(x))
    else:
        x_colors = fish_colors(len(x))
    ## y-axis
    y_max = kwargs.get('y_max', 'auto')
    if isinstance(y_max, numbers.Number):
        pass
    else:
        y_max = get_y_axis_max(x, chr, start=start, end=end)
    y_min = kwargs.get('y_min', 0)
    if y_min == 'auto' or isinstance(y_min, numbers.Number):
        pass
    else:
        y_min = 0
    ## bins
    nbins = kwargs.get('nbins', 700)
    ###########################################################################
    ## ini for region
    msg_list = []
    for bw,color in zip(x, x_colors):
        args_bw = {
            'color': color,
            'min_value': y_min,
            'max_value': y_max,
            'nbins': nbins,
        }
        msg = bw_to_config_f1(bw, **args_bw)
        if msg:
            msg_list.append(msg)
    ## add extra
    msg_list += ['[spacer]', '[x-axis]']
    config_text = '\n'.join(msg_list)
    ## save to file
    config_name = kwargs.get('config', 'config.ini')
    config_file = os.path.join(project_dir, config_name)
    with open(config_file, 'wt') as w:
        w.write(config_text+'\n')
    ###########################################################################
    ## run command
    plot_stdout = os.path.join(project_dir, 'stdout.log')
    plot_stderr = os.path.join(project_dir, 'stderr.log')
    plot_cmd = os.path.join(project_dir, 'run.sh')
    plot_file = kwargs.get('plot_file', None)
    if isinstance(plot_file, str):
        check_path(os.path.dirname(plot_file))
        pass
    else:
        plot_file = os.path.join(project_dir, 
                                 '{}_{}_{}.png'.format(chr, start, end))
    ## command
    cmd = ' '.join([
        '{}'.format(shutil.which('pyGenomeTracks')),
        '--tracks {}'.format(config_file),
        '--region {}:{}-{}'.format(chr, start, end),
        '-o {}'.format(plot_file),
        '--dpi 150',
        '--trackLabelFraction 0.2',
        '1> {}'.format(plot_stdout),
        '2> {}'.format(plot_stderr),
    ])
    with open(plot_cmd, 'wt') as w:
        w.write(cmd+'\n')
    if file_exists(plot_file):
        out = None
    else:
        out = cmd
    return out



def bw_to_tracks(x, **kwargs):
    """Generate the plots

    Parameters
    ----------
    x : list
        list of bigWig files
        
    Keyword parameters
    ------------------
    outdir: Directory to save the config.ini and plot
            if `None`, set the `cwd()`.
    region_list: regions in BED3 format
    y_min:  The minimum value on y-axis, default: 0
    y_max:  The maximum value on y-axis, default: 'auto'
    colors: Set colors for the tracks/bigWig files, 
            if `auto`, use colors by `fishualize` package
            from R, or set a list of colors by name, 
            `['red', 'green', 'blue']`, the number should
            match the bigWig files
    config: The name of `ini` file, defualt: `outdir/config.ini`
    plot_frefix: The name of the plot file, 'chr_start_end.png'
    """
    ## regions
    region_list = kwargs.get('region_list', None)
    if file_exists(region_list):
        regions = load_regions_from_bed(region_list)
    else:
        regions = load_regions_from_bw(x[0]) # chr:(start:end)
    ## run
    i = 0
    ia = len(regions)
    for chr,pos in regions.items():
        i += 1
        args_k = {
            'start': pos[0],
            'end': pos[1],
        }
        kwargs.update(args_k)
        cmd = bw_to_config(x, chr, **kwargs)
        print('{}/{} - {}'.format(i, ia, chr))
        if isinstance(cmd, str):
            os.system(cmd)



def te_track_plot(x, outdir=None):
    """Create plots for each TE
    
    Parameters
    ----------
    x : list
        list of project_dirs (pipe)
        
    outdir : str
        direcotory to save the config.ini and plots
        
    cmd:
    for single chr, requirements:
    1. config.ini # min/max/...
    2. cmd.sh
    
    pyGenomeTracks --tracks tracks.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o master_bigwig.png
    """
    # outdir
    outdir = file_abspath(outdir)
    # bam list
    bam_01 = get_te_bam(x[0]) # for headers
    bam_regions = load_regions_from_bam(bam_01)
    # bam to bw
    bw_list = [get_te_bw(i) for i in x] # list of list
    bw_list = [bw for sub in bw_list for bw in sub if bw] # remove None
    # config
    ini = os.path.join(outdir, 'config.ini')
    s = get_config(bw_list, ini)
    # command list
    cmd_list = []
    for k,v in bam_regions.items():
        plot_out = os.path.join(outdir, k+'.png')
        plot_stdout = os.path.join(outdir, 'stdout.log')
        plot_stderr = os.path.join(outdir, 'stderr.log')
        cmd = ' '.join([
            'pyGenomeTracks',
            '--tracks {}'.format(ini),
            '--region {}'.format(v),
            '-o {}'.format(plot_out),
            '--dpi 150',
            '--trackLabelFraction 0.2',
            '1> {}'.format(plot_stdout),
            '2> {}'.format(plot_stderr),
        ])
        if not file_exists(plot_out):
            cmd_list.append(cmd)
        # break
    # save command
    cmd_all = os.path.join(outdir, 'run_tracks.sh')
    with open(cmd_all, 'wt') as w:
        w.write('\n'.join(cmd_list)+'\n')
    # run
    run_shell_cmd(cmd_all)
    # os.system(cmd_all)
        

def te_track(x, outdir):
    """Generate te tracks
    
    Parameters
    ----------
    x : str
        a file, saving the (project_dirs), one per line 
        
    outdir : str
        directory, saving the plots
    """
    with open(x) as r:
        prj_list = r.readlines()
    prj_list = [i.strip() for i in prj_list]
    # output
    check_path(outdir)
    te_track_plot(prj_list, outdir)
    
        
        

def fish_colors(n=4, colorPal=1):
    """Pick colors

    Parameters
    ----------
    n : int
        number of colors return, default: 4
        
    colorPal : int or str
        choose the group of colors, 1 to N, or the name of fish
        candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
        'Centropyge_loricula'], default: [1]

    # colors from fishualize R package
    # https://nschiett.github.io/fishualize/articles/overview_colors.html

    Colors from R package: 
    > fishualize::fish_pal(option = "Scarus_hoefleri")(10)
     [1] "#D2372CFF" "#E25719FF" "#F0780BFF" "#F89814FF"
     [5] "#EFB52BFF" "#C6CB43FF" "#7AD45CFF" "#00CE7DFF"
     [9] "#00BCABFF" "#0499EAFF"

    > fishualize::fish_pal(option = "Gramma_loreto")(10)
     [1] "#020122FF" "#1E2085FF" "#4029CBFF"
     [4] "#6628EEFF" "#901CEDFF" "#B804CAFF"
     [7] "#D61693FF" "#E6445DFF" "#EE7A30FF"
    [10] "#F0BF0BFF"

    > fishualize::fish_pal(option = "Centropyge_loricula")(10)
     [1] "#8F1D1EFF" "#B30029FF" "#DF002AFF"
     [4] "#FF7D1AFF" "#FFBD17FF" "#E7BE5AFF"
     [7] "#988591FF" "#0043A0FF" "#001A72FF"
    [10] "#000000FF"
    """
    # pre-defined colors
    c1 = ['#D2372C', '#E25719', '#F0780B', '#F89814', '#EFB52B', 
        '#C6CB43', '#7AD45C', '#00CE7D', '#00BCAB', '#0499EA']
    c2 = ['#020122', '#1E2085', '#4029CB', '#6628EE', '#901CED', 
        '#B804CA', '#D61693', '#E6445D', '#EE7A30', '#F0BF0B']
    c3 = ['#8F1D1E', '#B30029', '#DF002A', '#FF7D1A', '#FFBD17', 
        '#E7BE5A', '#988591', '#0043A0', '#001A72', '#000000']
    # RGB (10 colors, repeat twice, 20 items)
    color_d = {
        'Scarus_hoefleri': c1*2,
        'Gramma_loreto': c2*2,
        'Centropyge_loricula': c3*2,
    }
    # get the fish_list
    fish_list = list(color_d.keys())
    # determine the fish (colorBy)
    if isinstance(colorPal, int):
        colorPal = colorPal - 1 # to 0-indexed
        if not colorPal in range(len(fish_list)):
            colorPal = 0
        fish = fish_list[colorPal]
    elif isinstance(colorPal, str):
        fish = colorPal if colorPal in color_d else 'Scarus_hoefleri'
    else:
        fish = 'Scarus_hoefleri'
    # output colors
    return color_d.get(fish, c1)[:n]



def bw_stats(x, chr=None, **kwargs):
    """Extract the max in the region, for each bam
    
    Parameters
    ----------
    x : str
        bigwig file
        
    chr : str 
        The chromosome name, if `None`, pick the first chromosome

    Keyword arguments:
            start: Starting position
            end:   Ending position
            type:  Summary type (mean, min, max, coverage, std), default 'mean'.
            nBins: Number of bins into which the range should be divided before
                   computing summary statistics. The default is 1.
            exact: By default, pyBigWig uses the same method as Kent's tools from UCSC
                   for computing statistics. This means that 'zoom levels' may be
                   used, rather than actual values (please see the pyBigWig repository
                   on github for further information on this). To avoid this behaviour,
                   simply specify 'exact=True'. Note that values returned will then
                   differ from what UCSC, IGV, and similar other tools will report.

    >>> import pyBigWig
    >>> bw = pyBigWig.open("test/test.bw")
    >>> bw.stats("1", 0, 3)
    [0.2000000054637591]
    pyBigWig.open(bw).stats('Idefix', 1, 7411, type='max')    
    """
    bw = pyBigWig.open(x)
    if chr is None: # choose chromosome
        chr = list(bw.chroms().keys()).pop(0)
    s = bw.stats(chr, **kwargs)
    bw.close()
    return s
    
    

def get_y_axis_max(x, chr=None, start=0, end=0, fix_end=True):
    """
    
    Parameters
    ----------
    x : list
        A list of bigwig files
        
    chr : str
        The name of chromosome 
        
    start : int 
        The start position 
        
    end : int
        The end position, 0 indicate the end of the chromosome
        
    fix_end : bool
        Fix the value, to nearest int (by breaks 5)
        
    """ 
    # file exists
    if isinstance(x, list):
        pass
    elif isinstance(x, str):
        x = [x]
    else:
        print('expect list, got {}'.format(type(x).__name__))
        return 'auto'
    x = [f for f in x if file_exists(f)]
    if len(x) == 0:
        print('no bigwig files found')
        return 'auto'
    # determine the args
    args = {
        'type': 'max',
    }
    if start > 0:
        args['start'] = start
    if end > 0:
        args['end'] = end
    # stats
    s = [bw_stats(bw, chr=chr, **args) for bw in x]
    smax = round(max(s).pop(0), 1)
    if fix_end:
        # breaks: by 5, 10
        breaks = 5
        smax = (int(round(smax/breaks, 1))+1)*breaks
    return smax
    
    
        
class TackPlot(object):
    """
    1. create track.ini
    2. CLI:
    $ pyGenomeTracks  --tracks tracks.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o master_bigwig.png
    """
    def __init__(self, **kwargs):
        pass
    
    
    def init_args(self):
        args_init = {
            'outdir': None
        }


    def get_template(self):
        """
        [track-01]
        file = 
        title = 
        height = 2
        color = 
        min_value = auto
        max_value = auto
        number_of_bins = 700
        file_type = bigwig
        [track-02]
        ...
        [spacer]
        [x-axis]
        """
        pass
    
    
    
    def get_colors(self):
        """
        Using R package: fishualize 
        """
        pass
    
    
    
    def run(self):
        pass
        
        
        
    def prep_ini():
        """Prepare for tracks.ini, 
        required arguments:
          - file:
          - color
          - type: [line|points|fill]
          - summary_method: [min|max]
          - min_value: [auto|numbers] 
          - max_value: [auto|numbers]
          - number_of_bins: 300
          - overlay_previous: [share-y|yes]
          - alpha: [0-1]
          - orientation: [inverted]


        optional arguments:
          - file_type: 
          - abc

        """
        pass


def main():
    if len(sys.argv) < 3:
        print('Usage: create_bw_tracks.py <prj.list> <outdir>')
        sys.exit(1)
    # args
    prj_txt = sys.argv[1]
    outdir = sys.argv[2]
    # load prj_list
    with open(prj_txt) as r:
        prj_list = r.readlines()
    prj_list = [i.strip() for i in prj_list]
    bw_list = [get_te_bw(i) for i in prj_list]
    bw_list = [i for sub in bw_list for i in sub if i]
    bw_to_tracks(bw_list, outdir=outdir)
    
    
if __name__ == '__main__':
    main()
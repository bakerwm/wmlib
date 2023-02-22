#!/usr/bin/env python

"""
Read Transposon from FlyBase (version FB2020_03, 2020-06-16)
url: ftp://ftp.flybase.org/releases/FB2020_03/precomputed_files/transposons/transposon_sequence_set.embl.txt.gz

Description
Input: EMBL
ouptut: GTF gene/CDS/LTR
output: BED LTR/non-LTR
output: fasta
output: txt, id conversion table

WARNING:
1. only extract Dmel records (126)
2. use 'dbxref' instead
3. features, manual curation required !!!! eg: P-element

# meta data
te_metadata_table.csv

embl_id: ID
fb_id: DR,ID
source_id: AC
flybase_name: DR,secondary
species: DR,secondary
te_name: DR,secondary
te_type_gff
te_type_fasta
te_subtype_fasta

# to-do
parse all features, saved as GTF
"""

import os
import sys
import re


class Embl(object):
    """
    Parse information from EMBL file
    multiple records
    """
    def __init__(self, input):
        self.input = input
        self.records = self.get_blocks()
        self.id = [i.id for i in self.get_blocks()] # DME9736
        self.fb_id = [i.fb_id for i in self.get_blocks()] # FBte0000104
        self.species = [i.species for i in self.get_blocks()] # Dmel
        self.fb_symbol = [i.fb_symbol for i in self.get_blocks()] # Idefix
        self.acc = [i.acc for i in self.get_blocks()] # AJ009736
        self.dbxref = [i.dbxref for i in self.get_blocks()] #
#         self.seq = self.get_seq()


    def get_blocks(self):
        """
        Read the EMBL record,
        split into records by "//" at the end
        """
        with open(self.input) as r:
            blocks = r.read().strip().split('\n//')
            for block in blocks:
                if len(block) > 0:
                    record = EmblR1(block)
                    # yield record # full records: 179, update: 2021-04
                    if record.species == 'Dmel':
                        yield record # Dmel records: 126, update: 2021-04


    def to_fa(self, output=None):
        if output is None:
            w = sys.stdout
        else:
            w = open(output, 'wt')
        for record in self.records:
            w.write('>{}\n{}'.format(record.fb_symbol, record.seq) + '\n')
        if not output is None:
            w.close()


    def to_gtf(self, output=None):
        if output is None:
            w = sys.stdout
        else:
            w = open(output, 'wt')
        for record in self.records:
            w.write('\n'.join(record.to_gtf()) + '\n')
        if not output is None:
            w.close()


    def to_gff(self, output=None):
        # determine output
        if output is None:
            w = sys.stdout
        else:
            w = open(output, 'wt')
        # iterator
        for record in self.records:
            w.write('\n'.join(record.to_gff()) + '\n')
        # close file
        if not output is None:
            w.close()


    def to_table(self, output=None):
        """
        Exract te_type
        output:
        embl_id, fb_id, source_id, flybase_name, species, te_name, te_type_gff,
        te_type_fasta, te_subtype_fasta
        """
        w = sys.stdout if output is None else open(output, 'wt')
        for record in self.records:
            w.write('\n'.join(record.to_table()) + '\n')
        if not output is None:
            w.close()


class EmblR1(object):
    """
    Extract the following info from single EMBL record (separated by "//")
    ID, name, flybase_id, flybase_name, source, CDS, LTR
    check if LTR or not: LTR, non_LTR
    x str a block of EMBL record
    begin with 'ID  '
    ID
    symbol
    name
    """
    def __init__(self, x):
        self.x = x # input
        self.info = self.read()
        self.id = self.get_id()
        self.fb_id = self.get_flybase('id')
        self.fb_symbol = self.get_flybase('symbol')
        self.species = self.get_flybase('species')
        self.acc = self.get_acc()
        self.dbxref = self.get_dbxref()
        self.seq = self.get_seq()
        self.features = self.get_features()
        self.size = len(self.seq)


    def valid_embl(self):
        """
        Check the input x, is valide EMBL
        criteria:
        1. begin with 'ID '
        2. start with '\W\W '
        """
        return self.x.startswith('ID   ') #


    def read(self):
        """
        split block by 'prefix'
        ID
        AC
        DR
        FT CDS/LTR/  source
        SQ
        """
        lines = self.x.strip().split('\n')
        # d save features
        ft_pre = None
        d = {} # key: ft, value = lines
        for line in lines:
            ft = line[:2] # feature type
            if re.search('^[A-Z]{2}', ft):
                ft_pre = ft # update pre
            else:
                if ft_pre:
                    ft = ft_pre # for seq lines
                else:
                    continue # skip blank lines
            d[ft] = d.get(ft, []) # init
            d[ft].append(line[5:]) # remove the prefix: "FT    "
        return d


    def get_flybase(self, group='symbol'):
        """
        arguments:
        group: str
            id, symbol, species in the line

        #!!!! extract flybase from "DR "
        example:
        DR   FLYBASE; FBte0000104; Dmel\Idefix.
        output: FBte0000104
        """
        lines = self.info.get('DR', None)
        if lines is None:
            out = None
        else:
            out = lines[0].rstrip('.') # first line
            fs = out.split('; ')
            if len(fs) == 3:
                fb, fb_id, fb_name = fs
                species, fb_symbol = fb_name.split('\\')
                if group == 'symbol':
                    out = fb_symbol
                elif group == 'id':
                    out = fb_id
                elif group == 'species':
                    out = species
                else:
                    out = None
        return out


    def get_id(self):
        """
        # for this file only !!!!
        example:
        ID   DME9736    standard; DNA; INV; 7411 BP.
        #!!! for this file only
        use dbxref FLYBASE as the name/id
        """
        lines = self.info.get('ID', None)
        if lines is None:
            id = None
        else:
            id = lines[0].split()[0]
        return id


    def get_acc(self):
        """
        d, a dict from read_block()
        example:
        AC   AJ009736;

        """
        lines = self.info.get('AC', None)
        if lines is None:
            acc = ''
        else:
            # only the first one for accession
            line = lines[0].rstrip(';')
            acc = line.split()[0]
        return acc


    def get_dbxref(self):
        """
        DR line         primary, secondary
        example:
        DR   FLYBASE; FBte0000104; Dmel\Idefix.
        """
        lines = self.info.get('DR', None) # get flybase only (db_xref)
        out = []
        if lines is None:
            pass
        elif isinstance(lines, list) and len(lines) > 0:
            for line in lines:
                line = line.rstrip('.').replace(';', '') # DR   FLYBASE FBte0000104 Dmel\\Idefix
                ts = line.split()
                if len(ts) >= 2:
                    out.append(':'.join(ts[:3]))
        return out


    def get_seq(self):
        """
        example:
        SQ   Sequence 7411 BP; 3047 A; 1363 C; 1109 G; 1892 T; 0 other;
             GTGACATATC CATAAGTCCC TAAGACTTAA GCATATGCCT ACATACTAAT ACACTTACAA        60
         """
        lines = self.info.get('SQ', None)
        if lines is None:
            seq = ''
        elif isinstance(lines, list) and len(lines) > 0:
            lines = lines[1:] # remove the first one
            lines = [re.sub('\s|[0-9]', '', i.rstrip()) for i in lines] # remove, white|digits
            seq = ''.join(lines)
        else:
            seq = ''
        return seq


    def get_features(self):
        """
        parse the features from read_features()
        example: [not a standard feature type]
        FT   source          AJ009736:1..7411
        FT   SO_feature      five_prime_LTR ; SO:0000425:1..600
        FT   SO_feature      three_prime_LTR ; SO:0000426:6841..7411
        FT   SO_feature      CDS ; SO:0000316:<988..2031
        """
        # return self.read_features()
        d = self.read_features() #
        out = []
        if isinstance(d, dict):
            for k, v in d.items():
                f = self.get_feature(v)
                # out.append(self.get_feature(v))
                if f:
                    out.append(f)
        return out


    def read_features(self):
        """
        split 'FT ... ' into features, saved in dict
        example: [not a standard feature type]
        FT   source          AJ009736:1..7411
        FT   SO_feature      five_prime_LTR ; SO:0000425:1..600
        FT   SO_feature      three_prime_LTR ; SO:0000426:6841..7411
        FT   SO_feature      CDS ; SO:0000316:<988..2031
        FT                   /name="Dmel\Idefix\gag"
        """
        lines = self.info.get('FT', None)
        if lines is None:
            return None
        # sub features
        d_sub = {}
        ft_pre = None
        for line in lines:
            # CDS ; SO:0000316:<469..1786
            # CDS ; SO:0000316:<1779..3772>
            line = line.replace(';', '')
            line = line.replace('<', '')
            line = line.replace('>', '')
            ts = line.rstrip().split(None, 2)
            if len(ts) == 1:
                # FT                   /name="Dmel\Idefix\pol"
                if ft_pre is None:
                    continue
                else:
                    ft = ft_pre
                    hit = ts[0]
            elif len(ts) == 2:
                # FT   source          AJ009736:1..7411
                key, ft = ts
                hit = ft # output
                ft_pre = ft
            elif len(ts) == 3:
                # FT   SO_feature      five_prime_LTR ; SO:0000425:1..600
                f, k, r = ts
                loc = r.split(':')[-1]
                ft = k + ':' + loc
                hit = ft
                ft_pre = ft
            else:
                continue
            d_sub[ft] = d_sub.get(ft, [])
            d_sub[ft].append(hit)
        return d_sub


    def get_feature(self, x, feature_list=[
        'exon', 'gene', 'CDS', 'mRNA', 'transcript',
        'five_prime_utr', 'three_prime_utr', 'start_codon',
        'five_prime_LTR', 'three_prime_LTR',
        'TATA_box']):
        """
        x string, block for each feature

        parse the name, db_xref, translation, ... #!!!

        FT   source          complement(AE002612:6924-6166)

        FT   SO_feature      three_prime_LTR ; SO:0000426:6841..7411
        FT   SO_feature      CDS ; SO:0000316:<988..2031
        FT                   /name="Dmel\Idefix\gag"
        FT                   /db_xref="FLYBASE:FBgn0027381"
        FT                   /db_xref="SPTREMBL:O96739"
        FT                   /db_xref="NCBI_PROTEIN:CAA08806.1"
        FT                   /translation="ARKLKDIMAVPQLSETHLNQLLNQIKELNYYDGAPGKLSGFVNQV
        output:
        (feature start end name db_xref translation)
        """
        ft = ''
        loc_start = ''
        loc_end = ''
        name = ''
        db_xref = ''
        translation = ''
        if isinstance(x, list) and len(x) > 0:
            pass
        else:
            return None
        # the feature_type, loc
        f = x[0]
        del x[0] # remove the first item
        p = re.search('(.*):<?(\d+)\.\.(\d+)>?', f)
        if p:
            ft, loc_start, loc_end = p.groups()
            if not ft in feature_list:
                return None
            if len(x) > 0:
                # get translation
                # name/db_xref/translation
                db_xref_list = []
                k_pre = None
                for s in x:
                    if s[0] == '/':
                        try:
                            k, v = s[1:].split('=', 1)
                        except:
                            print('!AAAA-2', s[0])
                        if k == 'name':
                            name = v
                        elif k == 'translation':
                            translation += v
                        elif k == 'db_xref':
                            db_xref_list.append(v)
                        else:
                            continue
                        k_pre = k
                    else:
                        if k_pre == 'translation':
                            translation += v
                        else:
                            continue
                # output
                db_xref = ','.join(db_xref_list)
        return [ft, loc_start, loc_end, name, db_xref, translation]


    def get_type_gff(self):
        """
        type in gff
        1. FT   SO_feature      terminal_inverted_repeat ; SO:0000481:1..27
        2. FT   SO_feature      non_LTR_retrotransposon ; SO:0000189 CC                   telomeric retrotransposon
        3. FT   SO_feature      five_prime_LTR ; SO:0000425:1..330
        4.
        """
        pass


    def get_type_fasta(self):
        """
        type in fasta
        """
        pass


    def get_subtype_fasta(self):
        """
        subtype in fasta
        """
        pass


    def to_gtf(self):
        """
        Convert to GTF format
        features: ft start end name db_xref,db_xref seq(aa)
        output:
            seqname
            source
            feature
            start
            end
            score
            strand
            frame
            group: gene_id ""; transcript_id ""; exon_number 1
        """
        gtf = []
        # gene
        attr1 = '; '.join([
            'gene_id \"{}\"'.format(self.fb_symbol),
            'gene_name \"{}\"'.format(self.fb_symbol),
            'gene_accession \"{}\"'.format(self.acc),
            'transcript_id \"{}\"'.format(self.id),
            'gene_source \"{}\"'.format('FlyBase'),
            'species \"{}\"'.format(self.species)
        ])
        gtf.append('\t'.join(
            [self.fb_symbol, 'FlyBase', 'gene', '1', str(self.size), '.', '+', '.', attr1]))
        # exon
        attr2 = '; '.join([
            'gene_id \"{}\"'.format(self.fb_symbol),
            'gene_name \"{}\"'.format(self.fb_symbol),
            'gene_accession \"{}\"'.format(self.acc),
            'transcript_id \"{}\"'.format(self.id),
            'exon_number \"{}\"'.format(1),
            'gene_source \"{}\"'.format('FlyBase'),
            'species \"{}\"'.format(self.species)
        ])
        gtf.append('\t'.join(
            [self.fb_symbol, 'FlyBase', 'exon', '1', str(self.size), '.', '+', '.', attr2]))
        # features
        for ft, start, end, name, dbxref, seq in self.features:
            if ft == '':
                continue # skip
            # name, "Dmel\Idefix\gag"
            name = name.replace('"', '').split('\\')[-1]
            if name == '':
                name = self.fb_symbol
            attr2 = '; '.join([
            'gene_id \"{}\"'.format(self.fb_symbol),
            'gene_name \"{}\"'.format(name),
            'transcript_id \"{}\"'.format(self.fb_symbol),
            'gene_source \"{}\"'.format('FlyBase')])
        return gtf


    def to_gff(self):
        """
        Convert to GFF format (GFF3)
        see: http://gmod.org/wiki/GFF3
        features: ft start end name db_xref,db_xref seq(aa)
        output:
            seqid
            source
            type
            start
            end
            score
            strand
            phase
            attributes
        for attributes:
            ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term
        example:
        ctg123 example mRNA            1050 9000 . + . ID=EDEN.1;Parent=EDEN;Name=EDEN.1;Index=1
        """
        gff = []
        # gene
        attr1 = ';'.join([
            'ID={}'.format(self.fb_symbol),
            'Name={}'.format(self.fb_symbol),
            'Accession={}'.format(self.acc),
            'Alias={}'.format(self.id),
            'Species={}'.format(self.species)
        ])
        gff.append('\t'.join(
            [self.fb_symbol, 'FlyBase', 'gene', '1', str(self.size), '.', '+', '.', attr1]))
        # features
        for ft, start, end, name, dbxref, seq in self.features:
            attr2 = ';'.join([
                'Parent=transcript',
                'ID={}'.format(self.symbol),
                'Name={}'.format(self.fb_symbol),
                'Dbxref={}'.format(dbxref),
                ])
            # gff.append('\t'.join(
            #     [self.name, 'FlyBase', ft, start, end, '.', '+', '.', attr2]))
        return gff


    def to_table(self):
        """
        # meta data
        te_metadata_table.csv
        embl_id: ID
        fb_id: DR,ID
        source_id: AC
        flybase_name: DR,secondary
        species: DR,secondary
        te_name: DR,secondary
        te_type_gff
        te_type_fasta
        te_subtype_fasta
        """
        # dbxref: DR: flybase: FBte..., Dmel/Idefix
        x = self.dbxref[0]
        db, fbid, fbname = x.split(':')
        species, te_name = fbname.split('\\')
        # output
        tab = '\t'.join([
            self.species,
            self.fb_id,
            self.fb_symbol,
            self.id,
            self.acc,
        ])
        return [tab]


def main():
    if len(sys.argv) < 3:
        sys.exit('read_embl.py [out_fmt:fa|gtf|gff|tab] <in.embl>')

    fmt = sys.argv[1].lower()
    input = sys.argv[2]
    output = sys.argv[3] if len(sys.argv) > 3 else None

    # read input
    p = Embl(input)

    if fmt == 'fa':
        p.to_fa(output)
    elif fmt == 'gff':
        p.to_gff(output)
    elif fmt == 'gtf':
        p.to_gtf(output)
    elif fmt == 'tab':
        p.to_table(output)
    else:
        sys.exit('unknown fmt: [fa|gtf|gff|tab], {}'.format(fmt))


if __name__ == '__main__':
    main()


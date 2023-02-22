#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download SRA data using aspera

## other tools
1. https://github.com/jhawkey/sra_read_downloader


example:

## from EBI
# command
$ ~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    -l 100M -T -k1 --mode recv \
    -P 33001 --host fasp.sra.ebi.ac.uk --user era-fasp \
    --file-list file.txt .
# notes
-l MAX-RATE, Max transfer rate, eg: 100M
-T Disable encryption 
-k Resume criterion: 0,3,2,1 
--mode MODE: send, recv 
-P SSH-PORT                     TCP port used for SSH authentication 
--host=HOSTNAME
--user=USERNAME
--file-list=FILENAME            File with list of sources
"""

import os
import sys
import pathlib
import argparse
import json
import requests
from hiseq.utils.utils import Config, run_shell_cmd, get_date, log
from hiseq.utils.file import file_exists, file_abspath, check_path, copy_file, fx_name

## simplifid version
## download by SRR id
##
## https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&fields=fastq_aspera,fastq_md5&format=json&query=run_accession="SRR11117696"

## step-1. query by srr
def query_srr_url(x):
    base_url = 'https://www.ebi.ac.uk/ena/portal/api/search?'
    fields=[
        'result=read_run',
        'format=json', 
        'fields=fastq_aspera,fastq_md5',
        'query=run_accession="{}"'.format(x)
    ]
    url = base_url + '&'.join(fields)
    return url


## step-2. validate url
def parse_url(url):
    response = requests.get(url)  # api call
    if response.ok:
        l = json.loads(response.text) # [{fastq_aspera,fastq_md5}, ...]
    else:
        l = []
    return l # list of dict
    

## step-3. download by aspera
def aspera_download(url, outdir=None, max_speed='100m'):
    """
    url: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/096/SRR11117696/SRR11117696_1.fastq.gz

    # example 
    ~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
        -l 100M -T -k1 --mode recv \
        -P 33001 url dest
    """
    if not isinstance(outdir):
        outdir = str(pathlib.Path.cwd())
    check_path(outdir)
    outfile = os.path.join(outdir, os.path.basename(url))
    if file_exists(outfile):
        log.info('file exists: {}'.format(outfile))
    else:
        uname = fx_name(url, fix_pe=False)
        cmd = ' '.join([
            '{}'.format(os.path.expanduser('~/.aspera/connect/bin/ascp')),
            '-i {}'.format(os.path.expanduer('~/.aspera/connect/etc/asperaweb_id_dsa.openssh')),
            '-l {}'.format(max_speed),
            '-T -k1 --mode recv -P 33001',
            'era-fasp@fasp.sra.ebi.ac.uk:{}'.format(url),
            outdir            
        ])
        cmd_txt = os.path.join(outdir, 'cmd'+uname+'.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # run
        try:
#             run_shell_cmd(cmd)
            print('download: {} files'.format(i))
        except:
            log.error('failed to run aspera ...')


def aspera_download2(url_list, outdir=None, max_speed='100m'):
    """
    url: fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/096/SRR11117696/SRR11117696_1.fastq.gz
    user: era-fasp

    # example 
    ~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
        -l 100M -T -k1 --mode recv \
        -P 33001 --host fasp.sra.ebi.ac.uk --user era-fasp \
        --file-list file.txt .
    """
    if not isinstance(outdir, str):
        outdir = str(pathlib.Path.cwd())
    outdir = file_abspath(outdir)
    check_path(outdir)
    # check exists files
    file_txt = os.path.join(outdir, 'file.txt')
    i = 0
    with open(file_txt, 'wt') as w:
        for f in url_list: 
            # url: fasp.sra.ebi.ac.uk:/vol1/fas...
            fname = fx_name(f, fix_pe=False)
            fout = os.path.join(outdir, fname)
            if not file_exists(fout):
                i += 1
                _, p = f.split(':')
                w.write(p+'\n')
    # prepare commands
    if i == 0:
        log.info('no urls in {}'.format(file_txt))
    else:
        stdout = os.path.join(outdir, 'aspera.stdout')
        stderr = os.path.join(outdir, 'aspera.stderr')
        cmd = ' '.join([
            '{}'.format('~/.aspera/connect/bin/ascp'),
            '-i {}'.format('~/.aspera/connect/etc/asperaweb_id_dsa.openssh'),
            '-l {}'.format(max_speed),
            '-T -k1 --mode recv -P 33001',
            '--host fasp.sra.ebi.ac.uk',
            '--user era-fasp',
            '--file-list {}'.format(file_txt),
            outdir,
            '1> {} 2> {}'.format(stdout, stderr)
        ])
        cmd_txt = os.path.join(outdir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # run
        try:
            run_shell_cmd(cmd)
#             print('download: {} files'.format(i))
        except:
            log.error('failed to run aspera ...')
 
        
## step-4: wrap all functions
def download_ena(x, outdir=None):    
    if not isinstance(outdir, str):
        outdir = str(pathlib.Path.cwd())
    outdir = file_abspath(outdir)
    check_path(outdir)
    # make sure, x is list
    if isinstance(x, str):
        x = [x]
    elif isinstance(x, list):
        pass
    else:
        raise ValueError('unknown x, expect str or list, got {}'.format(
            type(x).__name__))
    # filter by SRR
    x = [i for i in x if i.startswith('SRR')] # fixed, expand required
    f_list = [] # list of fastq_aspera files
    if len(x) > 0:
        # out, list
        url_list = [query_srr_url(i) for i in x]
        for u in url_list:
            l = parse_url(u) # list of dict {fastq_aspera,fastq_md5}
            for d in l:
                fq = d.get('fastq_aspera', None)
                if isinstance(fq, str):
                    if ';' in fq:
                        fq_list = fq.split(';')
                    else:
                        fq_list = [fq]
                    f_list.extend(fq_list)
    # download
    aspera_download2(f_list, outdir)
    

def main():
    # parse arguments
    if len(sys.argv) < 3:
        print('Usage: download_ena.py outdir [srr ...]')
        sys.exit(1)
    outdir = sys.argv[1]
    srr_list = sys.argv[2:]
    download_ena(srr_list, outdir)
        

        
        
if __name__ == '__main__':
    main()
    
# END

"""
ENA PortAPI

## Documentation
https://www.ebi.ac.uk/ena/portal/api/
https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access.html

## search from EBI (API)

curl -X GET --header 'Accept: application/json'
'https://www.ebi.ac.uk/ena/portal/api/search?<search_definition>'


https://www.ebi.ac.uk/ena/portal/api/search?result=read_run

query: &query=run_accession="SRR11117696"
get all available fields: "fields=all".
fields: &fields=all
output: &format=json


# Example
https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&fields=all&format=json&query=study_accession="PRJNA607400"

# Read files for {read_run} fileds
# return: json-array
# selected files

study_accession      : "PRJNA607400"
run_accession        : "SRR11117696"
submission_accession : "SRA1047238"
experiment_accession : "SRX7786259"
sample_accession     : "SAMN14168443"
sample_alias         : "GSM4333015"
experiment_alias     : "GSM4333015"
sample_title         : "GFP_Pol2_ChIP_rep2"
sample_description   : "GFP_Pol2_ChIP_rep2"
read_count           : "7108907"
base_count           : "2132672100"
fastq_aspera         : "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_2.fastq.gz"
fastq_ftp            : "ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_2.fastq.gz"
fastq_md5            : "fdbeeccd7677f8e91ccad45afe9215f4;76fa40bcc5e0bcfa7e137e649a34889a"
tax_id               : "10090"
scientific_name      : "Mus musculus"
instrument_platform  : "ILLUMINA"
instrument_model     : "HiSeq X Ten"
library_layout       : "PAIRED"
library_selection    : "ChIP"
library_strategy     : "ChIP-Seq"
library_source       : "GENOMIC"
study_title          : "Genome-wide analyses of chromatin interactions after the loss of Pol I, Pol II and Pol III"


# all available fileds

[
{"study_accession":"PRJNA607400","secondary_study_accession":"SRP250063","sample_accession":"SAMN14168443","secondary_sample_accession":"SRS6201758","experiment_accession":"SRX7786259","run_accession":"SRR11150179","submission_accession":"SRA1047238","tax_id":"10090","scientific_name":"Mus musculus","instrument_platform":"ILLUMINA","instrument_model":"HiSeq X Ten","library_name":"","library_layout":"PAIRED","nominal_length":"","library_strategy":"ChIP-Seq","library_source":"GENOMIC","library_selection":"ChIP","read_count":"7108907","base_count":"2132672100","center_name":"GEO","first_public":"2020-06-04","last_updated":"2020-06-04","experiment_title":"HiSeq X Ten sequencing; GSM4333015: GFP_Pol2_ChIP_rep2; Mus musculus; ChIP-Seq","study_title":"Genome-wide analyses of chromatin interactions after the loss of Pol I, Pol II and Pol III","study_alias":"PRJNA607400","experiment_alias":"GSM4333015","run_alias":"GSM4333015_r1","fastq_bytes":"641320931;570947734","fastq_md5":"fdbeeccd7677f8e91ccad45afe9215f4;76fa40bcc5e0bcfa7e137e649a34889a","fastq_ftp":"ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_2.fastq.gz","fastq_aspera":"fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_2.fastq.gz","fastq_galaxy":"ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/079/SRR11150179/SRR11150179_2.fastq.gz","submitted_bytes":"","submitted_md5":"","submitted_ftp":"","submitted_aspera":"","submitted_galaxy":"","submitted_format":"","sra_bytes":"","sra_md5":"","sra_ftp":"","sra_aspera":"","sra_galaxy":"","cram_index_ftp":"","cram_index_aspera":"","cram_index_galaxy":"","sample_alias":"GSM4333015","broker_name":"","nominal_sdev":"","first_created":"2020-06-04","sample_description":"GFP_Pol2_ChIP_rep2","parent_study":"PRJNA9559","accession":"SAMN14168443","bio_material":"","cell_line":"V6.5;mouse embryonic stem cells","cell_type":"","collected_by":"","collection_date":"","country":"","cultivar":"","culture_collection":"","description":"HiSeq X Ten sequencing; GSM4333015: GFP_Pol2_ChIP_rep2; Mus musculus; ChIP-Seq","dev_stage":"","ecotype":"","environmental_sample":"false","germline":"false","identified_by":"","isolate":"","isolation_source":"","location":"","mating_type":"","serotype":"","serovar":"","sex":"","submitted_sex":"","specimen_voucher":"","strain":"","sub_species":"","sub_strain":"","tissue_lib":"","tissue_type":"","variety":"","checklist":"","depth":"","elevation":"","altitude":"","environment_biome":"","environment_feature":"","environment_material":"","temperature":"","salinity":"","sampling_campaign":"","sampling_site":"","sampling_platform":"","protocol_label":"","project_name":"","host":"","host_tax_id":"","host_status":"","host_sex":"","submitted_host_sex":"","host_body_site":"","host_gravidity":"","host_phenotype":"","host_genotype":"","host_growth_conditions":"","environmental_package":"","investigation_type":"","experimental_factor":"","sample_collection":"","sequencing_method":"","target_gene":"","ph":"","sample_title":"GFP_Pol2_ChIP_rep2","sample_material":"","taxonomic_identity_marker":"","assembly_quality":"","assembly_software":"","taxonomic_classification":"","completeness_score":"","contamination_score":"","binning_software":"","lat":"","lon":"","sample_capture_status":"","collection_date_submitted":"","submission_tool":""}
]


## 2. Searchable fields

A full list of the latest result fields can be fetched from the following endpoint:
/searchFields?result=<resultId>

or by specific data portal:
/searchFields?result=<resultId>&dataPortal=pathogen

### 2.1 Avalilable <resultId> 

- sample fields:     [sample], [portal: ENA, ]
- read fields:       [read_run, read_experiment, read_study] 
- analysis fields:   [analysis, analysis_study]
- assembly fields:   [assembly]
- sequence fields:   [sequence_release, sequence_update]  
- Contig set fields: [wgs_set, tsa_set]  
- Coding fields:     [coding_release, coding_update] 
- Noncoding fields:  [noncoding_release, noncoding_update]


## 3. Returnable fields

A full list of the latest result fields can be fetched from the following endpoint:
/returnFields?result=<resultId>

or by specific data portal:
/returnFields?result=sample&dataPortal=pathogen


### 3.1 Available <resultId> 

- sample fields:     [sample], [portal: ENA, ]
- read fields:       [read_run, read_experiment, read_study] 
- analysis fields:   [analysis, analysis_study]
- assembly fields:   [assembly]
- sequence fields:   [sequence_release, sequence_update]  
- Contig set fields: [wgs_set, tsa_set]  
- Coding fields:     [coding_release, coding_update] 
- Noncoding fields:  [noncoding_release, noncoding_update]


## 4. Examples 

### 4.1 Search available fields

Fetch available results
Fetch the list of results that can be searched against in the pathogen data portal.
https://www.ebi.ac.uk/ena/portal/api/results?dataPortal=pathogen

Fetch searchable fields
Fetch the list of fields that can be searched for the assembly result.
https://www.ebi.ac.uk/ena/portal/api/searchFields?result=assembly

Fetch returnable fields
Fetch the list of fields that can be returned in the report for the analysis result
https://www.ebi.ac.uk/ena/portal/api/returnFields?result=analysis

### 4.2 Search for read data 

+ Search for read data using sample fields 

Find all public Salmonella read data collected in 2016. Return a list of the run accessions with FTP FASTQ file links. Include the MD5 checksums for the files.
https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=collection_date>=2016-01-01%20AND%20collection_date<=2016-12-31%20AND%20tax_tree(590)&fields=fastq_ftp,fastq_md5&limit=0

The query part of the url:
query=collection_date>=2016-01-01%20AND%20collection_date<=2016-12-31%20AND%20tax_tree(590)


https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=JSON&fields=fastq_aspera,fastq_md5&query=study_accession="PRJNA607400"

https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&fields=all&format=json&query=study_accession="PRJNA607400"



"""









"""
## Deprecated
anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra


## NCBI
$ cat file.txt
/refseq/release/viral/viral.1.1.genomic.fna.gz
/refseq/release/viral/viral.2.1.genomic.fna.gz

$ ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -T -k1 --mode recv --host ftp.ncbi.nlm.nih.gov --user anonftp --file-list file.lst .

## EBI
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -T -P33001 fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz .


ftp: ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/099/SRR11117699/SRR11117699_2.fastq.gz
asper: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR180/00


"""
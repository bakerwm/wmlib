#!/bin/bash 

# Check overlap between TE results

function ovl () {
    # BED6 format
    # a: query 
    # b: subject 
    # report number records a overlapped in b
    bedtools intersect -u -a $1 -b $2 | wc -l
}
export -f ovl

function fname () {
    # DNAseq_TE_001_F-element.notupd.merge.cutoff_3.bed -> F-element
    # gal4_ins.illumina.bed -> illumina
    # gal4_ins.ont.bed -> ont
    fn=$(basename $1)
    fn=${fn/.bed}
    fn=${fn/gal4_ins.}
    #
    fn=${fn/DNAseq_TE_}
    fn=${fn/.cutoff*}
    fn=${fn/.notupd.merge}
    echo ${fn}
}
export -f fname

function ovl2 () {
    # check overlapped records
    na=$(cat $1 | wc -l) # query
    nb=$(cat $2 | wc -l) # te
    nab=$(ovl $1 $2)
    nba=$(ovl $2 $1)
    #
    fa=$(fname $1)
    fb=$(fname $2)
    # extract only specific TE
    fte=${fa/*_}
    nte=$(grep -c "${fte}:" $2)
    # echo ${fte}
    echo ${fa},${fb},${na},${nb},${nte},${nab},${nba}
}
export -f ovl2

function main () {
    # scan cutoff
    # local ilmn=$(realpath -s data/te_ins/gal4_ins.illumina.bed)
    # local ont=$(realpath -s data/te_ins/gal4_ins.ont.bed)
    local ilmn="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.illumina.bed"
    local ont="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.ont.bed"
    # 1. ilmn
    ovl2 $1 ${ilmn}
    ovl2 $1 ${ont}
}
export -f main

# file header
echo fileA,fileB,numA,numB,numTE,AinB,BinA
parallel -j 1 main ::: results/*/*/*cutoff_*.bed

function ovl () {
    # BED6 format
    # a: query 
    # b: subject 
    # report number records a overlapped in b
    bedtools intersect -u -a $1 -b $2 | wc -l
}
export -f ovl

function fname () {
    # DNAseq_TE_001_F-element.nodup.merge.cutoff_3.bed -> F-element
    # gal4_ins.illumina.bed -> illumina
    # gal4_ins.ont.bed -> ont
    fn=$(basename $1)
    fn=${fn/.bed}
    fn=${fn/gal4_ins.}
    #
    fn=${fn/DNAseq_TE_}
    fn=${fn/.cutoff*}
    fn=${fn/.nodup.merge}
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

function ovl3 () {
    # scan cutoff
    local ilmn="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.illumina.bed"
    local ont="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.ont.bed"
    # 1: merged.bed
    # 2: cutoff
    input=$1
    cutoff=$2
    fc="${input/.bed}.cutoff_$2.bed"
    awk -v i=$2 '$5>i' ${input} > ${fc}
    s1=$(ovl2 ${fc} ${ilmn})
    s2=$(ovl2 ${fc} ${ont})
    echo ${cutoff},${s1}
    echo ${cutoff},${s2}
}
export -f ovl3

# scan cutoff
parallel -j 1 ovl3 {1} {2} ::: results/*2_roo*/*/*merge.bed ::: $(seq 10 2 100)
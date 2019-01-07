#!/usr/bin/env bash
# run STAR
#

[[ $# -lt 3 ]] && echo "Usage: map_star.sh <prj> <out_path> <r1>" && exit 0
star_idx="$HOME/software/piPipes/current/common/hg19/STARIndex/"
# star_idx="/home/wangming/software/piPipes/piPipes-1.5.0/common/dm3/STARIndex/"
# star_idx="/home/wangming/data/genome/dm3/STAR_index"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3_transposon"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3_genome_transposon"

# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3_genome_transposon"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/STAR_index/dm3_transposon"

# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/dm3.genome_transposon_STAR_index"
# star_idx="/home/wangming/data/genome/dm3/STAR_index"
# star_idx="/home/wangming/data/genome/dm3/canonical_transposons/dm3.genome+transposon.STAR_index"
# star_idx="/home/wangming/software/piPipes/piPipes-1.5.0/common/dm3/STARIndex/bbbbb_STAR_index/"
# star_idx="/home/wangming/data/genome/dm3/STAR_index/bbbbb_STAR_index/"
# star_idx="/home/wangming/data/genome/dm3/dm3.genome_transposon.STAR_index"
cpu=40
file_reader="-"

function is_gzip() {
    n=$(file -b $1)
    tag=0
    [[ ${n} = gzip* ]] && tag=1
    echo ${tag}
}


function run_star(){
    prjname=$1
    outdir=$2
    r1=$3
    # r2=$4
    STAR --runMode alignReads \
        --genomeDir ${star_idx} \
        --readFilesIn ${r1} \
        --readFilesCommand zcat \
        --runThreadN ${cpu} \
        --genomeLoad LoadAndRemove \
        --limitBAMsortRAM 10000000000 \
        --limitOutSAMoneReadBytes 1000000 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${outdir}/${prjname}.\
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --outFilterScoreMin 0 \
        --outFilterScoreMinOverLread 0.72 \
        --outFilterMatchNmin 0 \
        --outFilterMatchNminOverLread 0.72 \
        --outFilterMultimapScoreRange 1 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMax 0 \
        --alignIntronMin 21 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMunmapped None \
        --outReadsUnmapped Fastx \
        --outSJfilterReads Unique \
        --seedSearchStartLmax 20 \
        --seedSearchStartLmaxOverLread 1.0 \
        --chimSegmentMin 0
       
    # # output file: ${outdir}/${prjname}.Aligned.sortedByCoord.out.bam
    # samtools view -@ ${cpu} -bhS -F 4 ${outdir}/${prjname}.Aligned.out.sam | \
    #     samtools sort -n -@ ${cpu} -o ${outdir}/${prjname}.Aligned.out.bam -
    samtools index ${outdir}/${prjname}.Aligned.sortedByCoord.out.bam
}


# parameters
prjname=$1
outdir="$2/$prjname"
r1=$(readlink -fs $3) 

mkdir -p $outdir
run_star $prjname $outdir $r1

#

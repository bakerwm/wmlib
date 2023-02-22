# hiseq trim -1 data/raw_data/*1.fq.gz -2 data/raw_data/*2.fq.gz -m 20 --recursive --cut-after-trim 7 -p 12 -o data/clean_data -a AGATCGGAAGAGCACA -A CTGTCTCTTATACACA
# hiseq trim -1 data/raw_data/*1.fq.gz -2 data/raw_data/*2.fq.gz -m 20 --recursive --cut-after-trim 7 -p 12 -o data/clean_data_v2 -a AGATCGGAAGAGCACA -A CTGTCTCTTATACACA --cut-to-length 90
# hiseq qc -i data/clean_data_v2/*/*gz -o data/clean_data_v2/qc -p 8 -j 2 -g 
# hiseq qc -i data/clean_data/*/*gz -o data/clean_data/qc -p 4 -j 2 -g 
# hiseq qc -i data/raw_data/*gz -o data/raw_data/qc -p 4 -j 2 -g

# clean
# hiseq trim -1 data/raw_data/*1.fq.gz -2 data/raw_data/*2.fq.gz -m 20 --recursive --cut-after-trim 7,-7 --cut-to-length 97 -o data/clean_data -a AGATCGGAAGAGCACACGTCTGAAC -A CTGTCTCTTATACACATCTCCGAG -p 16 --times 4 -j 8

# hiseq align -a bowtie2 -p 20 -X "--no-mixed --no-discordant --very-fast-local" --index-list data/db/te/te -1 data/clean_data/*/*1.fq.gz -2 data/clean_data/*/*2.fq.gz -o align
# hiseq align -a bowtie2 -p 20 -X "--no-mixed --no-discordant --very-fast-local" --index-list data/db/dm6/dm6 -1 data/clean_data/*/*1.fq.gz -2 data/clean_data/*/*2.fq.gz -o align

#### overlapped with ONT/Illumina

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
    # local ilmn=$(realpath -s data/te_ins/gal4_ins.illumina.bed)
    # local ont=$(realpath -s data/te_ins/gal4_ins.ont.bed)
    local ilmn="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.illumina.bed"
    local ont="/data/yulab/wangming/work/yu_2023/projects/20230217_zp_TE_yy209/data/te_ins/gal4_ins.ont.bed"
    # 1. ilmn
    ovl2 $1 ${ilmn}
    ovl2 $1 ${ont}
}
export -f ovl3

####
function get_te_region () {
    # fastq name: DNAseq_TE_2_roo_1.fq.gz
    fq1=$1
    fn=$(basename ${fq1/_1.fq.gz})
    fn=${fn/DNAseq_TE*_} # TE name
    # te size
    ts=$(realpath data/db/te/te.fa.fai)
    s=$(awk -v i=${fn} '$1==i {print $2}' ${ts})
    echo ${fn}:1-${s}
}
export -f get_te_region

function find_ins () {
    local fq1=$1
    local out_dir=$2
    # local fq2=$2
    # local te_region=$4 # eg: roo:1-9092
    local fq2=${fq1/_1.fq/_2.fq}
    local te_region=$(get_te_region ${fq1})
    local prefix=$(basename ${fq1/_1.fq.gz})
    fq1=$(realpath ${fq1})
    fq2=$(realpath ${fq2})
    out_dir=$(realpath ${out_dir})
    echo ">>> Processing ${prefix}"

    # 1. Remove paired-mapped reads on TEs
    echo "1. Exclude paired-end reads on TE"
    index1=$(realpath "data/db/te/te")
    out_dir1="${out_dir}/${prefix}/1.paired_to_te"
    bam1="${out_dir1}/${prefix}.bam"
    log1="${out_dir1}/log.bowtie2.out"
    unmap1="${out_dir1}/${prefix}.unmap.fq"
    cmd1="${out_dir1}/cmd.sh"
    [[ ! -d ${out_dir1} ]] && mkdir -p ${out_dir1}
    txt="bowtie2 -x ${index1} -p 12 --very-fast-local \
        --no-discordant --no-mixed \
        --no-unal --un-conc ${unmap1} \
        -1 ${fq1} -2 ${fq2} 2> ${log1} | \
        samtools view -bhS -o ${bam1} -"
    echo ${txt} > ${cmd1}
    [[ ! -f ${bam1} ]] && bash ${cmd1}
    unfq1="${unmap1/.fq}.1.fq"
    unfq2="${unmap1/.fq}.2.fq"

    # 2. Map read2 to TE
    echo "2. Map read2 to TE"
    index2=$(realpath "data/db/te/te")
    out_dir2="${out_dir}/${prefix}/2.read2_to_te"
    bam2="${out_dir2}/${prefix}.bam"
    log2="${out_dir2}/log.bowtie2.out"
    unmap2="${out_dir2}/${prefix}.unmap.fq"
    cmd2="${out_dir2}/cmd.sh"
    [[ ! -d ${out_dir2} ]] && mkdir -p ${out_dir2}
    txt="bowtie2 -x ${index2} -p 12 --very-fast-local \
        --no-unal -U ${unfq2} 2> ${log2} | \
        samtools view -bhS - | \
        samtools sort -o ${bam2} -; \
        samtools index ${bam2}"
    echo ${txt} > ${cmd2}
    [[ ! -f ${bam2} ]] && bash ${cmd2}

    # 3. Exclude read2 not from specific TE
    echo "3. Exclude read2 not from specific TE"
    out_dir3="${out_dir}/${prefix}/3.read2_filter"
    bam3="${out_dir3}/${prefix}.bam"
    bam3_fq1="${out_dir3}/${prefix}.r1.fq"
    bam3_fq2="${out_dir3}/${prefix}.r2.fq"
    bam3_rid="${out_dir3}/${prefix}.r2.rid.txt"
    cmd3="${out_dir3}/cmd.sh"
    [[ ! -d ${out_dir3} ]] && mkdir -p ${out_dir3}
    # -F 0x4  exclude read2 unmapped
    # -f 0x10 include read2 on reverse strand of TE
    txt="samtools view -bhS -F 0x4 -f 0x10 -@ 8 ${bam2} ${te_region} | \
        samtools sort -@ 8 -o ${bam3} -  \n
        samtools index ${bam3} \n
        samtools fastq -0 ${bam3_fq2} ${bam3} \n
        samtools view ${bam3} | cut -f 1 > ${bam3_rid} \n
        seqkit grep -f ${bam3_rid} ${unfq1} > ${bam3_fq1}"
    echo -e ${txt} > ${cmd3}
    [[ ! -f ${bam3_fq1} ]] && bash ${cmd3}

    # 4. Map read1 to genome
    echo "4. Map read1 to genome"
    index4=$(realpath "data/db/dm6/dm6")
    out_dir4="${out_dir}/${prefix}/4.read1_to_genome"
    bam4="${out_dir4}/${prefix}.bam"
    log4="${out_dir4}/log.bowtie2.out"
    cmd4="${out_dir4}/cmd.sh"
    [[ ! -d ${out_dir4} ]] && mkdir -p ${out_dir4}
    # run
    txt="bowtie2 -x ${index4} -p 12 --very-fast-local \
        --no-unal -U ${bam3_fq1} 2> ${log4} | \
        samtools view -@ 8 -bhS - | \
        samtools sort -@ 8 -o ${bam4} - \n"
    echo -e ${txt} > ${cmd4}
    [[ ! -f ${bam4} ]] && bash ${cmd4}

    # 5. Filter read1
    echo "5. Filter read1 on genome"
    out_dir5="${out_dir}/${prefix}/5.filter_read1"
    bam5="${out_dir5}/${prefix}.bam"
    bam5_nodup="${out_dir5}/${prefix}.nodup.bam"
    bed5="${out_dir5}/${prefix}.nodup.bed"
    bed5_merge="${out_dir5}/${prefix}.nodup.merge.bed"
    bed5_filt="${out_dir5}/${prefix}.nodup.merge.cutoff_3.bed"
    log5="${out_dir5}/log.sambamba.out"
    cmd5="${out_dir5}/cmd.sh"
    [[ ! -d ${out_dir5} ]] && mkdir -p ${out_dir5}
    txt="# filter by mapQ \n
        samtools view -bhS -q 10 ${bam4} > ${bam5}
        # filt PCR duplicates \n
        sambamba markdup -r -t 8 ${bam5} ${bam5_nodup} 2> ${log5}\n
        # convert to bed \n
        bedtools bamtobed -i ${bam5_nodup} > ${bed5} \n
        # merge bed \n
        bedtools merge -d 200 -c 4,5,6 -o first,count,first -i ${bed5} > ${bed5_merge} \n
        # filter by cutoff=3 \n
        awk '\$5>3' ${bed5_merge} > ${bed5_filt}"
    echo -e ${txt} > ${cmd5}
    [[ ! -f ${bed5_filt} ]] && bash ${cmd5}

    # 6. Stat
    echo "6. Summary"
    stat6="${out_dir}/${prefix}/${prefix}.stat.csv"
    if [[ ! -f ${stat6} ]] ; then
        n_total=$(grep 'reads; of these' ${log1} | awk '{print $1}')
        n_te=$(samtools view -c -f 16 ${bam1})
        n_r2=$(cat ${bam3_rid} | wc -l)
        n_r1_genome=$(samtools view -c ${bam4})
        n_r1_q10=$(samtools view -c ${bam5})
        n_r1_nodup=$(samtools view -c ${bam5_nodup})
        n_r1_filt=$(awk '{s+=$5}END{print s}' ${bed5_filt})
        # stat
        n_r2_unmap=$((${n_total} - ${n_te} - ${n_r2}))
        n_r1_unmap=$((${n_r2} - ${n_r1_genome}))
        n_r1_multi=$((${n_r1_genome} - ${n_r1_q10}))
        n_r1_dup=$((${n_r1_q10} - ${n_r1_nodup}))
        n_r1_cutoff=$((${n_r1_nodup} - ${n_r1_filt}))
        # for insertions
        n_ins_merge=$(cat ${bed5_merge} | wc -l)
        n_ins_filt=$(cat ${bed5_filt} | wc -l)
        n_ins_cutoff=$((${n_ins_merge} - ${n_ins_filt}))
        # report
        echo name,total,te_paired,r2_unmap,r1_unmap,r1_multi,r1_dup,r1_cutoff,r1_output,ins_merge,ins_cutoff,ins_output > ${stat6}
        echo ${prefix},${n_total},${n_te},${n_r2_unmap},${n_r1_unmap},${n_r1_multi},${n_r1_dup},${n_r1_cutoff},${n_r1_filt},${n_ins_merge},${n_ins_cutoff},${n_ins_filt} >> ${stat6}
    fi

    # 7. Overlapped
    echo "7. Overlapped with ONT,Illumina TE insertions"
    stat7="${out_dir}/${prefix}/${prefix}.overlap.csv"
    if [[ ! -f ${stat7} ]] ; then
        echo fileA,fileB,numA,numB,numTE,AinB,BinA > ${stat7}
        ovl3 ${bed5_filt} >> ${stat7}
    fi

    # Finish
    echo ""
}
export -f find_ins



parallel -j 4 find_ins {} results ::: data/clean_data/*001*/*1.fq.gz





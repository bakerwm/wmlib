#!/usr/bin/env bash

# Name: run_para.sh
# Date: 2020-01-15 
# Author: MW


[[ $# -lt 1 ]] && echo "Usage: $(basename $0) <cmds.txt> <num>" && exit 0
cmd_txt=$(readlink -fs $1)
num=$2
[[ -z ${num} ]] && num=4 #
re='^[0-9]+$' # only integer allowed
[[ ${num} =~ $re ]] || (echo "${n} - non-integer found, exiting!" && exit 1)

## run in parallel
function run_parallel(){
    plist=()
    cat $1 | while read line ; do
        echo $line | bash & # run 
        plist+=("$!")
    done
    wait ${plist[@]}
}

## split cmd file in lines
## -l num, number of lines 
function split_txt(){
    txt_in=$1
    num=$2
    split -a 30 -l $num $txt_in cmd_split_ # split
}

## switch to directory of cmd_txt
cmd_dir=$(dirname ${cmd_txt})
pushd ${cmd_dir}

## remove old files: cmd_split_aa*
old_check=$(ls -1 | grep -c ^cmd_split_aaaaaaaaaa)
[[ ${old_check} -gt 0 ]] && rm cmd_split_aaaaaaaaaa* 

## split cmd
cmd_txt_name=$(basename ${cmd_txt})
split_txt ${cmd_txt_name} ${num}

## run all commands
for cmd_sub in cmd_split_aa* ; do
    run_parallel ${cmd_sub} ${num}
    echo $cmd_sub && sleep 3 # between sub commands, waiting for 30s
done

## rm sub files
rm cmd_split_aaaaaaaaaa*

## go back
popd 


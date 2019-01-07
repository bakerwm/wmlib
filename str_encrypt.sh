#!/bin/bash 

# encrypt and decrypt strings
#
# example
# bash str_encrypt.sh <string> 0|1
#
#


function usage() {
    cat <<EOF
Usage: $0 [option] <string>

Options:
    -r    decrypt string

EOF
}

function encrypt() {
    echo "$1" | base64 -i
}

function decrypt() {
    echo "$1" | base64 -d
}

decrypt=0
while getopts "hr:" OPTION ; do
    case $OPTION in 
        h)    usage && exit 0 ;;
        r)    decrypt=1 ;;
    esac
done
str=$1


if [[ $decrypt -gt 0 ]] ; then
    decrypt $str
else
    encrypt $str
fi

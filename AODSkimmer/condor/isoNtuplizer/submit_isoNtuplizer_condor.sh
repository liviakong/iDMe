#!/bin/bash

flist=$1
nev=$2

if [ "$#" -ne 2 ]; then
    echo "Bad arguments: need to provide fileList(.txt) and nEvents"
    exit
fi

flist_full=`realpath $flist`
fname=`echo $flist_full | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1`

condor_submit isoNtuplizer_config.jdl -append "Arguments = ${fname} ${nev}" -append "transfer_input_files = ${flist_full}"
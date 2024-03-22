#!/bin/bash

flist=$1
year=$2

flist_full=`realpath $flist`
echo ${flist_full}
fname=`echo $flist_full | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1`

mkdir -p Logs

condor_submit filterEffGetter_config.jdl -append "Arguments = ${fname} ${year}" -append "transfer_input_files = ${flist_full}"
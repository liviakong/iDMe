#!/bin/bash

flist=$1
year=$2
nThreads=$3
isData=$4
isSignal=$5
topDirName=$6

flist_full=`realpath $flist`
fname=`echo $flist_full | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1`

outDirName="${fname}"

mkdir -p split_fileLists
mkdir -p Logs

cp ${flist_full} .
split -d -l 20 --additional-suffix ".txt" ${fname}.txt ${fname}_
rm ${fname}.txt
for sublist in `ls ${fname}_*.txt`
do
	mv $sublist split_fileLists/$sublist
	sublist_name=`echo $sublist | cut -d "." -f 1`
	sublist_full=`realpath split_fileLists/$sublist`
	condor_submit ntuplizer_config_generic.jdl -append "Arguments = ${sublist_name} ${year} ${nThreads} ${isData} ${isSignal} ${outDirName} ${topDirName}" -append "transfer_input_files = ${sublist_full}" -append "request_cpus = ${nThreads}"
done

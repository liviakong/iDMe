#!/bin/bash

flist=$1
config=$2
outCollName=$3
nThreads=$4
config=`realpath $config`
flist_full=`realpath $flist`
fname=`echo $flist_full | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1`
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
	condor_submit reMiniAOD_config.jdl -append "Arguments = ${sublist_name} ${outCollName} ${nThreads}" -append "transfer_input_files = ${sublist_full},${config}"
done

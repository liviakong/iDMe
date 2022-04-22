#!/bin/bash
eval `scram runtime -sh`
dir=$1
name=`echo $dir | rev | cut -d "/" -f 1 | rev`
xrd="root://cmsxrootd.fnal.gov/"
touch files.txt
for f in `eos root://cmseos.fnal.gov/ ls $dir/*AOD*.root`
do
    echo $xrd$dir/$f >> files.txt
done
nfiles=`cat files.txt | wc -l`
if [ "$nfiles" -gt 50 ]
then
	split -d -l 50 --additional-suffix ".txt" files.txt merged_${name}_
	rm files.txt
	for sublist in `ls merged_${name}_*.txt`
	do
        out_fname=`echo $sublist | cut -d "." -f 1`
		edmCopyPickMerge inputFiles_load=$sublist outputFile="file:${out_fname}.root"
        xrdcp ${out_fname}.root "root://cmseos.fnal.gov/${dir}/${out_fname}.root"
        rm $sublist
        rm ${out_fname}.root
	done
else
	edmCopyPickMerge inputFiles_load=files.txt outputFile="file:merged_${name}.root"
    xrdcp merged_${name}.root "root://cmseos.fnal.gov/${dir}/merged_${name}.root"
fi


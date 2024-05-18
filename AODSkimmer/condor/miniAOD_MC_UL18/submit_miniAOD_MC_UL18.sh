#!/bin/bash

flist=$1
nThreads=$2
nsplit=$3

flist_full=`realpath $flist`
fname=`echo $flist_full | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1`

mass=`echo $fname | cut -d "_" -f 1-2`
ctau=`echo $fname | cut -d "_" -f 3`
outDirName="${mass}/${ctau}"

mkdir -p split_fileLists
mkdir -p Logs
outPath=/store/group/lpcmetx/iDMe//Samples/signal/2018/MINIAOD/${outDirName}/
xrdfs root://cmseos.fnal.gov/ mkdir -p ${outPath}

cp ${flist_full} .
split -d -l ${nsplit} --additional-suffix ".txt" ${fname}.txt ${fname}_
rm ${fname}.txt
for sublist in `ls ${fname}_*.txt`
do
	mv $sublist split_fileLists/$sublist
	sublist_name=`echo $sublist | cut -d "." -f 1`
	sublist_full=`realpath split_fileLists/$sublist`
	condor_submit miniAOD_MC_UL18_config.jdl -append "Arguments = ${sublist_name} ${nThreads} ${outPath}" -append "transfer_input_files = ${sublist_full},/uscms_data/d3/sbrightt/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/AODSkimmer/scripts/MC_miniAOD_UL18_condor.py" -append "request_cpus = ${nThreads}"
done

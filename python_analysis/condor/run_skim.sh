#!/bin/bash
jobname_base=$1
jobname=$2
mode=$3
n_cores=$4
outDir=$5

tar xzf ${jobname}.tar.gz
python condor_skim.py $jobname_base $jobname $mode $n_cores $outDir
#echo "xrdcp -f *.root root://cmseos.fnal.gov/${outDir}/"
#xrdcp -f *.root root://cmseos.fnal.gov/${outDir}/
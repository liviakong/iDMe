#!/bin/bash
jobname_base=$1
jobname=$2
mode=$3
n_cores=$4
outDir=$5
MET_cut=$6

source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh
tar xzf ${jobname}.tar.gz
python condor_skim_rdf.py $jobname_base $jobname $mode $n_cores $outDir $MET_cut
#echo "xrdcp -f *.root root://cmseos.fnal.gov/${outDir}/"
#xrdcp -f *.root root://cmseos.fnal.gov/${outDir}/
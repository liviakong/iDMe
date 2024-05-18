#!/bin/bash

fname=$1
nThreads=$2
outPath=$3

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
if ! [ -r CMSSW_10_6_28/src ] ; then
    scram p CMSSW CMSSW_10_6_28
fi

mv ${fname}.txt CMSSW_10_6_28/src/
mv MC_miniAOD_UL18_condor.py CMSSW_10_6_28/src/
cd CMSSW_10_6_28/src
eval `scram runtime -sh`
scram b -j 4
cmsRun MC_miniAOD_UL18_condor.py flist=${fname}.txt numThreads=${nThreads}
mv test_output.root MiniAOD_${fname}.root
xrdcp -f MiniAOD_${fname}.root root://cmseos.fnal.gov/${outPath}/MiniAOD_${fname}.root
echo "Copied MiniAOD_${fname}.root"
echo "Done"

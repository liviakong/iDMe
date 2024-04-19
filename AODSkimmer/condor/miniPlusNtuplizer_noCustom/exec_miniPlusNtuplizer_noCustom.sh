#!/bin/bash

fname=$1
year=$2
nThreads=$3
isData=$4
isSignal=$5
outPath=$6

xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26_noCustomMini.tar.gz .
tar -xzf ntuplizer_CMSSW_10_6_26_noCustomMini.tar.gz

mv ${fname}.txt CMSSW_10_6_26/src/iDMe/AODSkimmer
cd CMSSW_10_6_26/src/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram b ProjectRename
eval `scram runtime -sh`
cd iDMe/AODSkimmer
cmsRun miniPlusElectronNtuplizer_cfg.py flist=${fname}.txt data=${isData} signal=${isSignal} year=${year} numThreads=${nThreads}
mv test_output.root ntuples_${fname}.root
xrdcp -f ntuples_${fname}.root root://cmseos.fnal.gov/${outPath}/ntuples_${fname}.root
echo "Copied ntuples_${fname}.root"
echo "Done"

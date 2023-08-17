#!/bin/bash

fname=$1
nev=$2

if [ "$#" -ne 2 ]; then
    echo "Bad arguements: Need sampleName and nEvents"
    exit
fi

xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26.tar.gz .
tar -xzf ntuplizer_CMSSW_10_6_26.tar.gz

mv ${fname}.txt CMSSW_10_6_26/src/iDMe/AODSkimmer
cd CMSSW_10_6_26/src/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram b ProjectRename
eval `scram runtime -sh`
cd iDMe/AODSkimmer
cmsRun run_isoNtuplizer_cfg.py flist=${fname}.txt nEvents=${nev}
mv test_output.root isoNtuples_${fname}.root
xrdcp isoNtuples_${fname}.root root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//Samples/isoNtuples/
echo "Done"

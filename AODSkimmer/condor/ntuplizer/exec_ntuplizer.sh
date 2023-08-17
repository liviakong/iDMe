#!/bin/bash

fname=$1

xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26.tar.gz .
tar -xzf ntuplizer_CMSSW_10_6_26.tar.gz

mv ${fname}.txt CMSSW_10_6_26/src/iDMe/AODSkimmer
cd CMSSW_10_6_26/src/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram b ProjectRename
eval `scram runtime -sh`
cd iDMe/AODSkimmer
cmsRun run_AODntuplizer_cfg.py flist=${fname}.txt
mv test_output.root ntuples_${fname}.root
xrdcp ntuples_${fname}.root root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//Samples/Ntuples/
echo "Done"

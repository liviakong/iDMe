#!/bin/bash

fname=$1
year=$2

xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26.tar.gz .
tar -xzf ntuplizer_CMSSW_10_6_26.tar.gz

mv ${fname}.txt CMSSW_10_6_26/src/iDMe/AODSkimmer
cd CMSSW_10_6_26/src/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram b ProjectRename
eval `scram runtime -sh`
cd iDMe/AODSkimmer
cmsRun genFilterEfficiency_cfg.py flist=${fname}.txt
mv filter_eff.root filterEff_${fname}.root
xrdcp -f filterEff_${fname}.root root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//Samples/filterEffs/signal/${year}/filterEff_${fname}.root
echo "Copied filterEff_${fname}.root"
echo "Done"
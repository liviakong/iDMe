#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

BASE=/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/iDMeSkimmer
HEAD=/uscms/home/sbrightt/nobackup/iDMe/compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26

cd $HEAD

if ! [ -r CMSSW_10_6_26/src ] ; then
    scram p CMSSW CMSSW_10_6_26
fi

cd CMSSW_10_6_26/src
eval `scram runtime -sh`
rm -rf *
mkdir iDMeAnalysis
mkdir iDMeAnalysis/iDMeSkimmer
cp -r $BASE/plugins iDMeAnalysis/iDMeSkimmer/
cp -r $BASE/python iDMeAnalysis/iDMeSkimmer/
cp -r $BASE/test iDMeAnalysis/iDMeSkimmer/
cp $BASE/scripts/run_ntuplizer_cfg.py iDMeAnalysis/iDMeSkimmer/
cp -r /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/RecoVertex .

scram b -j 8

cd $HEAD
tar -czf ntuplizer_CMSSW_10_6_26.tar.gz CMSSW_10_6_26/
xrdcp -f ntuplizer_CMSSW_10_6_26.tar.gz root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/

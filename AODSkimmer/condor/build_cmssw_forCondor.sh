#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

BASE=/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/
HEAD=/uscms/home/sbrightt/nobackup/iDMe/compiled_CMSSW_envs/ntuplizer_CMSSW_10_6_26

mkdir -p $HEAD
cd $HEAD

if ! [ -r CMSSW_10_6_26/src ] ; then
    scram p CMSSW CMSSW_10_6_26
fi

cd CMSSW_10_6_26/src
eval `scram runtime -sh`
rm -rf *

mkdir iDMeAnalysis
mkdir iDMeAnalysis/iDMeSkimmer
mkdir iDMeAnalysis/AODSkimmer
mkdir iDMeAnalysis/CustomTools

cp -r $BASE/iDMeSkimmer/plugins iDMeAnalysis/iDMeSkimmer/
cp -r $BASE/iDMeSkimmer/python iDMeAnalysis/iDMeSkimmer/
cp -r $BASE/iDMeSkimmer/test iDMeAnalysis/iDMeSkimmer/

cp -r $BASE/AODSkimmer/plugins iDMeAnalysis/AODSkimmer/
cp -r $BASE/AODSkimmer/python iDMeAnalysis/AODSkimmer/
cp -r $BASE/AODSkimmer/test iDMeAnalysis/AODSkimmer/

cp -r $BASE/CustomTools/* iDMeAnalysis/CustomTools/

cp $BASE/iDMeSkimmer/scripts/run_ntuplizer_cfg.py iDMeAnalysis/iDMeSkimmer/
cp $BASE/AODSkimmer/scripts/run_AODntuplizer_cfg.py iDMeAnalysis/AODSkimmer/
cp $BASE/AODSkimmer/scripts/run_isoNtuplizer_cfg.py iDMeAnalysis/AODSkimmer/
cp $BASE/AODSkimmer/scripts/miniPlusElectronNtuplizer_cfg.py iDMeAnalysis/AODSkimmer/

cp -r /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/RecoVertex .
cp -r /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/RecoEgamma .
cp -r /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/PhysicsTools .

scram b -j 8

cd $HEAD
tar -czf ntuplizer_CMSSW_10_6_26.tar.gz CMSSW_10_6_26/
xrdcp -f ntuplizer_CMSSW_10_6_26.tar.gz root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/
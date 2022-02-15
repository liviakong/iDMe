#!/bin/bash

export BASEDIR=`pwd`

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820

tar xzf submit.tar.gz
cd CMSSW_10_6_26/src
scram b ProjectRename -j 4
eval `scram runtime -sh`

cd $BASEDIR

xrdcp -r CMSSW_10_6_26/src/ root://cmsxrootd.fnal.gov//store/group/lpcmetx/iDMe/test/

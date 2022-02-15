#! /bin/bash

## This script is used to re-RECO/re-AOD existing DIGIRAWHLT signal samples 
## with the ECAL RecHit timing cleaning turned OFF. Study to see if removing 
## timing cleaning improves displaced electron reconstruction efficiency

## Usage: ./runOffGridpack.sh DIGIRAW_FILE (with EOS path)

export BASEDIR=`pwd`
FILE_IN=$1
FNAME_IN=`echo ${FILE_IN} | rev | cut -d "/" -f 1 | rev`
nevent=1000
echo "Processing file ${FNAME_IN}"

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820
xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/cmssw_10_6_26_noTimingClean.tar.gz .
tar xzf cmssw_10_6_26_noTimingClean.tar.gz
cd CMSSW_10_6_26/src
scram b ProjectRename
eval `scram runtime -sh`

xrdcp root://cmsxrootd.fnal.gov/$FILE_IN .

# AOD production
echo "1.) Generating AOD"
cmsDriver.py step2 \
    --filein file:${FNAME_IN} \
    --fileout file:AOD_${FNAME_IN} \
    --mc --eventcontent AODSIM --datatier AODSIM --runUnscheduled \
    --conditions 106X_upgrade2018_realistic_v16_L1v1 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI \
    --procModifiers premix_stage2 \
    --nThreads 1 --era Run2_2018 --python_filename AOD_cfg.py --no_exec \
    --customise Configuration/DataProcessing/Utils.addMonitoring -n ${nevent} || exit $?;
cmsRun -p AOD_cfg.py

# MINIAOD production
echo "2.) Generating MINIAOD"
cmsDriver.py step3 \
        --filein file:AOD_${FNAME_IN} \
        --fileout file:MINIAOD_${FNAME_IN} \
        --mc --eventcontent MINIAODSIM --datatier MINIAODSIM --runUnscheduled \
        --conditions 106X_upgrade2018_realistic_v16_L1v1 --step PAT \
        --procModifiers run2_miniAOD_UL --nThreads 8 --era Run2_2018 \
        --python_filename MINIAOD_cfg.py --no_exec \
        --geometry DB:Extended --customise Configuration/DataProcessing/Utils.addMonitoring \
        -n ${nevent} || exit $?;
cmsRun -p MINIAOD_cfg.py

echo "Done"

#! /bin/bash

## This script is used to produce AOD files from a gridpack for
## 2018 data. The CMSSW version is 10_2_3 and all four lifetimes are
## produced: {1, 10, 100, 1000} mm per seed.
##
## The lifetime replacement no longer occurs at the LHE level (i.e.
## manually replacing the lifetime in LHE events) but rather at the
## Pythia hadronizer level. For private production, there are four 
## different hadronizers, one for each ctau, which gets called by this
## script as appropriate. For official central production, can have the
## calling script `sed` into the one hadronizer file to change the 
## lifetime accordingly.
##
## Currently MINIAOD production is commented out to save time (and we don't use it).

## Usage: ./runOffGridpack.sh gridpack_file.tar.xz

export BASEDIR=`pwd`
GP_f=$1
GRIDPACKDIR=${BASEDIR}/gridpacks
LHEDIR=${BASEDIR}/mylhes
SAMPLEDIR=${BASEDIR}/samples
[ -d ${LHEDIR} ] || mkdir ${LHEDIR}

HADRONIZER="iDMe_externalLHEProducer_and_PYTHIA8_Hadronizer"
namebase=${GP_f/.tar.xz/}
nevent=1000
#nevent=10

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820
if ! [ -r CMSSW_10_2_16_UL/src ] ; then
  scram p CMSSW CMSSW_10_2_16_UL
fi
if ! [ -r CMSSW_10_6_26/src ] ; then
    scram p CMSSW CMSSW_10_6_26
fi
cd CMSSW_10_6_26/src
eval `scram runtime -sh`
scram b -j 4
tar xaf ${GRIDPACKDIR}/${GP_f}
sed -i 's/exit 0//g' runcmsgrid.sh
ls -lrth

RANDOMSEED=`od -vAn -N4 -tu4 < /dev/urandom`
#Sometimes the RANDOMSEED is too long for madgraph
RANDOMSEED=`echo $RANDOMSEED | rev | cut -c 3- | rev`

echo "0.) Generating LHE"
sh runcmsgrid.sh ${nevent} ${RANDOMSEED} 4
namebase=${namebase}_$RANDOMSEED
cp cmsgrid_final.lhe ${LHEDIR}/${namebase}.lhe
echo "${LHEDIR}/${namebase}.lhe" 
rm -rf *
cd ${BASEDIR}

export SCRAM_ARCH=slc7_amd64_gcc820
if ! [ -r CMSSW_10_6_26/src ] ; then
    scram p CMSSW CMSSW_10_6_26
fi
cd CMSSW_10_6_26/src
rm -rf *
mkdir -p Configuration/GenProduction/python/

for ctau_mm in 1 10 100 1000
#for ctau_mm in 1
do
    cp "${BASEDIR}/conf/${HADRONIZER}_ctau-${ctau_mm}.py" Configuration/GenProduction/python/
    eval `scram runtime -sh`
    scram b -j 4
    echo "1.) Generating GEN-SIM for lifetime ${ctau_mm}"
    genfragment=${namebase}_GENSIM_cfg_ctau-${ctau_mm}.py
    cmsDriver.py Configuration/GenProduction/python/${HADRONIZER}_ctau-${ctau_mm}.py \
        --filein file:${LHEDIR}/${namebase}.lhe \
        --fileout file:${namebase}_GENSIM_ctau-${ctau_mm}.root \
        --mc --eventcontent RAWSIM --datatier GEN-SIM \
        --conditions 106X_upgrade2018_realistic_v16_L1v1 --beamspot Realistic25ns13TeVEarly2018Collision \
        --step GEN,SIM --era Run2_2018 --nThreads 1 \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --python_filename ${genfragment} --no_exec -n ${nevent} || exit $?;

    #Make each file unique to make later publication possible
    linenumber=`grep -n 'process.source' ${genfragment} | awk '{print $1}'`
    linenumber=${linenumber%:*}
    total_linenumber=`cat ${genfragment} | wc -l`
    bottom_linenumber=$((total_linenumber - $linenumber ))
    tail -n $bottom_linenumber ${genfragment} > tail.py
    head -n $linenumber ${genfragment} > head.py
    echo "    firstRun = cms.untracked.uint32(1)," >> head.py
    echo "    firstLuminosityBlock = cms.untracked.uint32($RANDOMSEED)," >> head.py
    cat tail.py >> head.py
    mv head.py ${genfragment}
    rm -rf tail.py

    cmsRun -p ${genfragment}

    # Step1 is pre-computed, since it takes a while to load all pileup pre-mixed samples
    # So just replace template with correct filenames and number of events
    echo "2.) Generating DIGI-RAW for lifetime ${ctau_mm}"
    cp ${BASEDIR}/DIGIRAW_template_2018_UL.py .
    sed -i "s/file:placeholder_in.root/file:${namebase}_GENSIM_ctau-${ctau_mm}.root/g" DIGIRAW_template_2018_UL.py
    sed -i "s/file:placeholder_out.root/file:${namebase}_DIGIRAW_ctau-${ctau_mm}.root/g" DIGIRAW_template_2018_UL.py
    sed -i "s/input = cms.untracked.int32(10)/input = cms.untracked.int32(${nevent})/g" DIGIRAW_template_2018_UL.py
    mv DIGIRAW_template_2018_UL.py ${namebase}_DIGIRAW_cfg_ctau-${ctau_mm}.py

    #echo "2.) Generating DIGI-RAW-HLT for lifetime ${ctau_mm}"
    #cmsDriver.py step1 \
    #    --filein file:${namebase}_GENSIM_ctau-${ctau_mm}.root \
    #    --fileout file:${namebase}_DIGIRAWHLT_ctau-${ctau_mm}.root \
    #    --era Run2_2018 --conditions 106X_upgrade2018_realistic_v16_L1v1 \
    #    --mc --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:@relval2018 \
    #    --procModifiers premix_stage2 \
    #    --datamix PreMix \
    #    --datatier GEN-SIM-DIGI-RAW --eventcontent PREMIXRAW \
    #    --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer17PrePremix-PUAutumn18_102X_upgrade2018_realistic_v15-v1/GEN-SIM-DIGI-RAW" \
    #    --number ${nevent} \
    #    --geometry DB:Extended --nThreads 1 \
    #    --python_filename ${namebase}_DIGIRAWHLT_cfg_ctau-${ctau_mm}.py \
    #    --customise Configuration/DataProcessing/Utils.addMonitoring \
    #    --no_exec || exit $?;
    cmsRun -p ${namebase}_DIGIRAW_cfg_ctau-${ctau_mm}.py

    #Go to older CMSSW for HLT
    cd ${BASEDIR}/CMSSW_10_2_16_UL/src
    eval `scram runtime -sh`
    scram b -j 4
    echo "3.) Generating HLT for lifetime ${ctau_mm}"
    cmsDriver.py step1 \
        --filein file:${BASEDIR}/CMSSW_10_6_26/src/${namebase}_DIGIRAW_ctau-${ctau_mm}.root \
        --fileout file:${namebase}_DIGIRAWHLT_ctau-${ctau_mm}.root \
        --python_filename ${namebase}_DIGIRAWHLT_cfg_ctau-${ctau_mm}.py \
        --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier GEN-SIM-RAW --conditions 102X_upgrade2018_realistic_v15 \
        --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
        --step HLT:2018v32 --geometry DB:Extended --era Run2_2018 --no_exec --mc -n ${nevent} || exit $?;
    cmsRun -p ${namebase}_DIGIRAWHLT_cfg_ctau-${ctau_mm}.py
    cp ${namebase}_DIGIRAWHLT_ctau-${ctau_mm}.root ${BASEDIR}/CMSSW_10_6_26/src

    #Go back to new CMSSW for AOD
    export SCRAM_ARCH=slc7_amd64_gcc820
    cd ${BASEDIR}/CMSSW_10_6_26/src
    eval `scram runtime -sh`
    scram b -j 4
    echo "4.) Generating AOD for lifetime ${ctau_mm}"
    cmsDriver.py step2 \
        --filein file:${namebase}_DIGIRAWHLT_ctau-${ctau_mm}.root \
        --fileout file:${namebase}_AOD_ctau-${ctau_mm}_year-2018.root \
        --mc --eventcontent AODSIM --datatier AODSIM --runUnscheduled \
        --conditions 106X_upgrade2018_realistic_v16_L1v1 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI \
        --procModifiers premix_stage2 \
        --nThreads 1 --era Run2_2018 --python_filename ${namebase}_AOD_cfg_ctau-${ctau_mm}.py --no_exec \
        --customise Configuration/DataProcessing/Utils.addMonitoring -n ${nevent} || exit $?;
    cmsRun -p ${namebase}_AOD_cfg_ctau-${ctau_mm}.py

    # MINIAOD production is commented out
    echo "5.) Generating MINIAOD"
    cmsDriver.py step3 \
            --filein file:${namebase}_AOD_ctau-${ctau_mm}_year-2018.root \
            --fileout file:${namebase}_MINIAOD_ctau-${ctau_mm}_year-2018.root \
            --mc --eventcontent MINIAODSIM --datatier MINIAODSIM --runUnscheduled \
            --conditions 106X_upgrade2018_realistic_v16_L1v1 --step PAT \
            --procModifiers run2_miniAOD_UL --nThreads 8 --era Run2_2018 \
            --python_filename ${namebase}_MINIAOD_cfg_ctau-${ctau_mm}.py --no_exec \
            --geometry DB:Extended --customise Configuration/DataProcessing/Utils.addMonitoring \
            -n ${nevent} || exit $?;
    cmsRun -p ${namebase}_MINIAOD_cfg_ctau-${ctau_mm}.py

    pwd
    cmd="ls -arlth *.root"
    echo $cmd && eval $cmd

    echo "DONE."
done
echo "ALL Done"

#! /bin/bash

## This script is used to produce AOD files from a gridpack for
## 2018 data
##
## Based on the MC&I UL MC instructions found at https://exo-mc-and-i.gitbook.io/exo-mc-and-interpretation/how-to-sample-production-private
##
## Usage: ./gridpackGen_UL18.sh gridpack_file.tar.xz

export BASEDIR=`pwd`
GP_f=$1
nevent=$2
ctau=$3
nthreads=$4
year=2018

if [ "$#" -ne 4 ]; then
    echo "Wrong number of arguments!"
    echo "Usage: ./genFromGridpack_UL.sh gridpack nevents ctau nThreads"
    exit 1
fi

echo "Generating signal MC for gridpack ${GP_f}"
echo "Lifetime set to ${ctau} mm"

GRIDPACKDIR=${BASEDIR}
LHEDIR=${BASEDIR}/mylhes
[ -d ${LHEDIR} ] || mkdir ${LHEDIR}

HADRONIZER="iDMe_pythiaGenFragment.py"
namebase=${GP_f/.tar.xz/}

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820
if ! [ -r CMSSW_10_2_16_UL/src ] ; then
  scram p CMSSW CMSSW_10_2_16_UL
fi
if ! [ -r CMSSW_10_6_28/src ] ; then
    scram p CMSSW CMSSW_10_6_28
fi

cd CMSSW_10_6_28/src
eval `scram runtime -sh`
scram b -j 4
tar xaf ${GRIDPACKDIR}/${GP_f}
sed -i 's/exit 0//g' runcmsgrid.sh
ls -lrth

RANDOMSEED=`od -vAn -N4 -tu4 < /dev/urandom`
#Sometimes the RANDOMSEED is too long for madgraph
RANDOMSEED=`echo $RANDOMSEED | rev | cut -c 3- | rev`

echo "0.) LHE Step"
sh runcmsgrid.sh ${nevent} ${RANDOMSEED} ${nthreads}
namebase=${namebase}_$RANDOMSEED
cp cmsgrid_final.lhe ${LHEDIR}/${namebase}.lhe
echo "${LHEDIR}/${namebase}.lhe" 
rm -rf *
cd ${BASEDIR}

export SCRAM_ARCH=slc7_amd64_gcc820
if ! [ -r CMSSW_10_6_28/src ] ; then
    scram p CMSSW CMSSW_10_6_28
fi
cd CMSSW_10_6_28/src
rm -rf *
mkdir -p Configuration/GenProduction/python/

cp "${BASEDIR}/${HADRONIZER}" Configuration/GenProduction/python/
eval `scram runtime -sh`
scram b -j 4

# Running GEN  Step
echo "########################################################################################"
echo "########################################################################################"
echo "1.) GEN Step"
genfragment=${namebase}_GEN_cfg_ctau-${ctau}.py
cmsDriver.py Configuration/GenProduction/python/${HADRONIZER} \
    --filein file:${LHEDIR}/${namebase}.lhe \
    --fileout file:${namebase}_GEN_ctau-${ctau}_year-${year}.root \
    --mc --eventcontent RAWSIM --datatier GEN \
    --step GEN --geometry DB:Extended \
    --conditions 106X_upgrade2018_realistic_v4 --beamspot Realistic25ns13TeVEarly2018Collision \
    --era Run2_2018 --nThreads $nthreads \
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

# Running SIM Step
echo "########################################################################################"
echo "########################################################################################"
echo "2.) SIM Step"
genfragment=${namebase}_SIM_cfg_ctau-${ctau}.py
cmsDriver.py step1 \
    --filein file:${namebase}_GEN_ctau-${ctau}_year-${year}.root \
    --fileout file:${namebase}_SIM_ctau-${ctau}_year-${year}.root \
    --mc --eventcontent RAWSIM --datatier GEN-SIM \
    --step SIM --geometry DB:Extended \
    --conditions 106X_upgrade2018_realistic_v11_L1v1 --beamspot Realistic25ns13TeVEarly2018Collision \
    --era Run2_2018 --nThreads $nthreads \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --python_filename ${genfragment} --no_exec --runUnscheduled -n ${nevent} || exit $?;

cmsRun -p ${genfragment}

# Running DigiPremix Step
echo "########################################################################################"
echo "########################################################################################"
echo "3.) DIGIPremix Step"
genfragment=${namebase}_DIGIPremix_cfg_ctau-${ctau}.py
cp $BASEDIR/template_DIGIPremix_cfg_UL2018.py ./${genfragment}
sed -i -e "s/PLACEHOLDER_IN.root/${namebase}_SIM_ctau-${ctau}_year-${year}.root/g" ${genfragment}
sed -i -e "s/PLACEHOLDER_OUT.root/${namebase}_DIGIPremix_ctau-${ctau}_year-${year}.root/g" ${genfragment}
sed -i -e "s/NEVT/${nevent}/g" ${genfragment}
sed -i -e "s/NTHREAD/${nthreads}/g" ${genfragment}

#cmsDriver.py step1 \
#    --filein file:${namebase}_SIM_ctau-${ctau}_year-${year}.root \
#    --fileout file:${namebase}_DIGIPremix_ctau-${ctau}_year-${year}.root \
#    --mc --eventcontent PREMIXRAW --datatier GEN-SIM-DIGI \
#    --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --datamix PreMix \
#    --geometry DB:Extended \
#    --pileup_input dbs:/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL18_106X_upgrade2018_realistic_v11_L1v1-v2/PREMIX \
#    --conditions 106X_upgrade2018_realistic_v11_L1v1 \
#    --era Run2_2018 --nThreads $nthreads \
#    --customise Configuration/DataProcessing/Utils.addMonitoring \
#    --python_filename ${genfragment} --no_exec --runUnscheduled -n ${nevent} || exit $?;

cmsRun -p ${genfragment}

# Doing HLT step in CMSSW 10_2_16_UL
cd ${BASEDIR}/CMSSW_10_2_16_UL/src
eval `scram runtime -sh`
scram b -j 4
echo "########################################################################################"
echo "########################################################################################"
echo "4.) HLT Step"
genfragment=${namebase}_HLT_cfg_ctau-${ctau}.py
cmsDriver.py step1 \
    --filein file:${BASEDIR}/CMSSW_10_6_28/src/${namebase}_DIGIPremix_ctau-${ctau}_year-${year}.root \
    --fileout file:${namebase}_HLT_ctau-${ctau}_year-${year}.root \
    --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW \
    --step HLT:2018v32 --geometry DB:Extended \
    --conditions 102X_upgrade2018_realistic_v15 \
    --era Run2_2018 --nThreads $nthreads \
    --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --python_filename ${genfragment} --no_exec -n ${nevent} || exit $?;

cmsRun -p ${genfragment}
cp ${namebase}_HLT_ctau-${ctau}_year-${year}.root ${BASEDIR}/CMSSW_10_6_28/src

# Doing RECO/AOD step
cd ${BASEDIR}/CMSSW_10_6_28/src
eval `scram runtime -sh`
scram b -j 4
echo "########################################################################################"
echo "########################################################################################"
echo "5.) RECO/AOD Step"
genfragment=${namebase}_AOD_cfg_ctau-${ctau}.py
cmsDriver.py step1 \
    --filein file:${namebase}_HLT_ctau-${ctau}_year-${year}.root \
    --fileout file:${namebase}_AOD_ctau-${ctau}_year-${year}.root \
    --mc --eventcontent AODSIM --datatier AODSIM \
    --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --geometry DB:Extended \
    --conditions 106X_upgrade2018_realistic_v11_L1v1 \
    --era Run2_2018 --nThreads $nthreads \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --python_filename ${genfragment} --no_exec --runUnscheduled -n ${nevent} || exit $?;

cmsRun -p ${genfragment}

# Doing MINIAOD Step
#echo "########################################################################################"
#echo "########################################################################################"
#echo "6.) MINIAOD Step"
#genfragment=${namebase}_MINIAOD_cfg_ctau-${ctau}.py
#cmsDriver.py step1 \
#    --filein file:${namebase}_AOD_ctau-${ctau}_year-${year}.root \
#    --fileout file:${namebase}_MINIAOD_ctau-${ctau}_year-${year}.root \
#    --mc --eventcontent MINIAODSIM --datatier MINIAODSIM \
#    --step PAT --geometry DB:Extended \
#    --conditions 106X_upgrade2018_realistic_v16_L1v1 \
#    --era Run2_2018 --nThreads $nthreads \
#    --procModifiers run2_miniAOD_UL \
#    --python_filename ${genfragment} --no_exec --runUnscheduled -n ${nevent} || exit $?;

#cmsRun -p ${genfragment}

# Doing MINIAOD Step with run2_miniAOD_devel modifier to (maybe) get the energy-regressed low-pT electrons
#echo "########################################################################################"
#echo "########################################################################################"
#echo "7.) MINIAOD Step with run2_miniAOD_devel modifier"
#genfragment=${namebase}_MINIAOD_devel_cfg_ctau-${ctau}.py
#cmsDriver.py step1 \
#    --filein file:${namebase}_AOD_ctau-${ctau}_year-${year}.root \
#    --fileout file:${namebase}_MINIAOD_devel_ctau-${ctau}_year-${year}.root \
#    --mc --eventcontent MINIAODSIM --datatier MINIAODSIM \
#    --step PAT --geometry DB:Extended \
#    --conditions 106X_upgrade2018_realistic_v16_L1v1 \
#    --era Run2_2018 --nThreads $nthreads \
#    --python_filename ${genfragment} --no_exec --runUnscheduled -n ${nevent} || exit $?;

#cmsRun -p ${genfragment}

pwd
cmd="ls -arlth *.root"
echo $cmd && eval $cmd

remoteDIR="/store/group/lpcmetx/iDMe//Samples/signal/2018"
#xrdcp -vf ${namebase}_HLT_ctau-${ctau}_year-${year}.root root://cmseos.fnal.gov/$remoteDIR/DIGIRAWHLT/${namebase}_HLT_ctau-${ctau}_year-${year}.root
xrdcp -vf ${namebase}_AOD_ctau-${ctau}_year-${year}.root root://cmseos.fnal.gov/$remoteDIR/AOD/${namebase}_AOD_ctau-${ctau}_year-${year}.root
#xrdcp -vf ${namebase}_MINIAOD_ctau-${ctau}_year-${year}.root root://cmseos.fnal.gov/$remoteDIR/MINIAOD/${namebase}_MINIAOD_ctau-${ctau}_year-${year}.root
#xrdcp -vf ${namebase}_MINIAOD_devel_ctau-${ctau}_year-${year}.root root://cmseos.fnal.gov/$remoteDIR/MINIAOD_devel/${namebase}_MINIAOD_devel_ctau-${ctau}_year-${year}.root

echo "DONE."

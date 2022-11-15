#!/bin/bash

fname=$1
outCollName=$2
nThreads=$3

xrdcp root://cmseos.fnal.gov//store/group/lpcmetx/iDMe//compiled_CMSSW_envs/reMini_lowptIDRegressNoCut.tar.gz .
tar -xzf reMini_lowptIDRegressNoCut.tar.gz

mv ${fname}.txt CMSSW_10_6_28/src/
mv mini_config.py CMSSW_10_6_28/src/
cd CMSSW_10_6_28/src/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
scram b ProjectRename
eval `scram runtime -sh`
for f in `cat ${fname}.txt`
do
    cp mini_config.py miniRun.py
    sed -i -e "s|FILE_PLACEHOLDER|$f|g" miniRun.py
    sed -i -e "s|NTHREAD|${nThreads}|g" miniRun.py
    cmsRun miniRun.py
    outdir=`echo $f | sed -e "s|root://cmsxrootd.fnal.gov/||g" | rev | cut -d "/" -f 2- | rev | sed -e "s/AOD/${outCollName}/g"`
    outfname=`echo $f | rev | cut -d "/" -f 1 | rev | sed -e "s/AOD/${outCollName}/g"`
    mv output.root ${outfname}
    xrdfs root://cmseos.fnal.gov/ mkdir -p ${outdir}
    xrdcp ${outfname} root://cmseos.fnal.gov/${outdir}/${outfname}
done
echo "Done"

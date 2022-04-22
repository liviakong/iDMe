#!/bin/bash

prefix=$PWD

for cat in `ls ${prefix}/fileLists/MINIAOD`
do
    cmsRun ${prefix}/scripts/run_ntuplizer_cfg.py flist=${prefix}/fileLists/MINIAOD/${list}
    fname=`echo ${list} | cut -d "." -f 1`
    mv test_output.root ntuples_MINIAOD_${fname}.root 
done

#!/bin/bash

mchi=$1
dmchi=$2

prefix="/store/group/lpcmetx/iDMe/Samples/"
cmseos="root://cmseos.fnal.gov/"

if [ "$#" -ne 2 ]
then
    echo "Wrong arguments! Usage: ./eos_dirproducer mchi dmchi"
    exit 0
fi

for ct in 1 10 100 1000
do
    eos $cmseos mkdir -p $prefix/DIGIRAWHLT/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}
    eos $cmseos mkdir -p $prefix/AOD/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}
    eos $cmseos mkdir -p $prefix/MINIAOD/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}
done

for ct in 1 10 100 1000
do
    for f in `eos $cmseos ls $prefix/*Mchi-${mchi}_dMchi-${dmchi}*_DIGIRAWHLT_ctau-${ct}.root`
    do
        eos $cmseos mv $prefix/$f $prefix/DIGIRAWHLT/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}/$f
    done

    for f in `eos $cmseos ls $prefix/*Mchi-${mchi}_dMchi-${dmchi}*_AOD_ctau-${ct}_*.root`
    do
        eos $cmseos mv $prefix/$f $prefix/AOD/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}/$f
    done

    for f in `eos $cmseos ls $prefix/*Mchi-${mchi}_dMchi-${dmchi}*_MINIAOD_ctau-${ct}_*.root`
    do
        eos $cmseos mv $prefix/$f $prefix/MINIAOD/Mchi-${mchi}_dMchi-${dmchi}_ctau-${ct}/$f
    done
done
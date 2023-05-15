#!/bin/bash

basedir=$1

cmseos="root://cmseos.fnal.gov/"
xrd="root://cmsxrootd.fnal.gov/"

for d in `eos $cmseos ls $basedir`
do
    filename="${d}.txt"
    touch $filename
    for f in `eos $cmseos ls $basedir/$d/*.root`
    do
        echo $xrd$basedir/$d/$f >> $filename
    done
done
#!/bin/bash

for f in `ls *UL*.txt`;
do
    campaign=`echo $f | cut -d "." -f 1 | rev | cut -d "_" -f 1 | rev`;
    for ff in `cat $f`;
    do
        dName=`echo $ff | cut -d "/" -f 2`;
        fullName=${campaign}_${dName};
        dasgoclient --limit=100 --query="file dataset=${ff}" >> fileLists/${fullName}.txt;
    done
done

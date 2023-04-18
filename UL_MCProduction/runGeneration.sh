#!/bin/bash

gridpack=$1
ctau=$2
year=$3
nevents=$4
nthreads=$5

export HOME=$PWD

tar -xzf submit.tar.gz
sed -i -e "s/SET_CTAU/${ctau}/g" iDMe_pythiaGenFragment.py
sh genFromGridpack_UL${year}.sh ${gridpack} ${nevents} ${ctau} ${nthreads}
#cd $HOME
#rm -rf *

exit 0

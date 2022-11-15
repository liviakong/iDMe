#!/bin/bash

gridpack=$1
ctau=$2
year=$3
nevents=$4
nthreads=$5

export HOME=$PWDs

tar -xzf submit.tar.gz
sed -i -e "s/SET_CTAU/${ctau}/g" iDMmu_pythiaGenFragment_noFilters.py
sh genFromGridpack_muon_UL${year}.sh ${gridpack} ${nevents} ${ctau} ${nthreads}
cd $HOME
rm -rf *

exit 0

#!/bin/bash

gridpack=$1
ctau=$2
ctau_mm=$3
year=$4
nevents=$5
nthreads=$6

export HOME=$PWD

tar -xzf submit.tar.gz
sed -i -e "s/SET_CTAU/${ctau}/g" iDMe_pythiaGenFragment.py

echo "gridpack = ${gridpack}"
echo "ctau = ${ctau}"
echo "ctau_mm = ${ctau_mm}"
echo "year = ${year}"
echo "nevents = ${nevents}"
echo "nthreads = ${nthreads}"

sh genFromGridpack_UL${year}.sh ${gridpack} ${nevents} ${ctau_mm} ${nthreads}
#cd $HOME
#rm -rf *

exit 0

#!/bin/bash
eosdir=$1

for f in `eos root://cmseos.fnal.gov/ find ${eosdir} -name "*.root" -type f`
do
	if [[ $f == *.root ]] 
	then 
		file=`echo $f | cut -d "/" -f 4-`
		echo "/${file}"
	fi
done

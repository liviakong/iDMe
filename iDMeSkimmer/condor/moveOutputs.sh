#!/bin/bash
pref="/store/group/lpcmetx/iDMe/"
for f in `eos root://cmseos.fnal.gov/ ls  $pref/Samples/Ntuples/ntuples_Mchi*.root`
do
	mass=`echo $f | cut -d "_" -f 2-3`
	ct=`echo $f | cut -d "_" -f 4`
	eos root://cmseos.fnal.gov/ rm $pref/Samples/Ntuples/signal/2018/${mass}/${ct}/$f
	eos root://cmseos.fnal.gov/ mv $pref/Samples/Ntuples/$f $pref/Samples/Ntuples/signal/2018/${mass}/${ct}/
done

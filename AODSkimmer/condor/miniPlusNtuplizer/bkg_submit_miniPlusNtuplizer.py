#!/usr/bin/env python
import sys
import os
import json

flist = sys.argv[1]
year = sys.argv[2]
nThreads = sys.argv[3]
isData = sys.argv[4]
isSignal = sys.argv[5]
suffix = sys.argv[6]
if len(sys.argv) > 7:
    nev_target = int(sys.argv[7])
else:
    nev_target = 10000

fname = flist.split("/")[-1].split(".")[0]
sample = fname.split("_")[0]
subsample = "_".join(fname.split("_")[1:])
outDirName=f"{sample}/{subsample}/"

# make output directory
os.system(f"xrdfs root://cmseos.fnal.gov/ mkdir -p /store/group/lpcmetx/iDMe//Samples/Ntuples/bkg_condor_{suffix}/{year}/{outDirName}/")
outPath = f"/store/group/lpcmetx/iDMe//Samples/Ntuples/bkg_condor_{suffix}/{year}/{outDirName}/"

with open(flist,'r') as fin:
    files = fin.read().splitlines()

nev_count = 0
subLists = [[]]
for line in files:
    f = line.split(" ")[0]
    nev = int(line.split(" ")[1])
    subLists[-1].append(f)
    nev_count += nev
    if nev_count >= nev_target:
        subLists.append([])
        nev_count = 0
if subLists[-1] == []:
    waste = subLists.pop(-1)

os.system('mkdir -p split_fileLists')
os.system('mkdir -p Logs')

for i,s in enumerate(subLists):
    sublist_name = f"{sample}_{subsample}_{i}"
    with open(f"split_fileLists/{sublist_name}.txt","w") as fout:
        fout.write("\n".join(s))
    sublist_full = f"{os.getcwd()}/split_fileLists/{sublist_name}.txt"
    condor_cmd = f"condor_submit miniPlusNtuplizer_config.jdl -append \"Arguments = {sublist_name} {year} {nThreads} {isData} {isSignal} {outPath}\" -append \"transfer_input_files = {sublist_full}\" -append \"request_cpus = {nThreads}\""
    os.system(condor_cmd)

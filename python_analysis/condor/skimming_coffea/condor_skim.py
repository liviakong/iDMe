#!/usr/bin/env python
import coffea.util as util
import os, stat
import json
import sys
import multiprocessing as mp
from multiprocessing import Process, Value, Array, Manager
import numpy as np
import subprocess

def skim_files(fileList,sampConfig):
    for infile in fileList:
        skimmer = fileSkimmer(infile,sampConfig,"cuts.py",mode=sampConfig['type'])
        skimmer.skim()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Bad input!")
        print("Usage: ./condor_skim jobname_base jobname mode n_cores")

    num_cpus = mp.cpu_count()
    jobname_base = sys.argv[1]
    jobname = sys.argv[2]
    mode = sys.argv[3]
    n_cores = int(sys.argv[4])
    outDir = sys.argv[5]
    from analysisTools import fileSkimmer
    with open('samples.json','r') as fin:
        samps = json.load(fin)
    files = samps['fileset']
    files = [f for f in files if f not in samps['blacklist']]
    """sub_lists = [list(l) for l in np.array_split(files,num_cpus)]
    processes = []
    for sublist in sub_lists:
        processes.append(Process(target=skim_files,args=(sublist,samps)))
    for p in processes:
        p.start()
    for p in processes:
        p.join()"""
    for inFile in files:
        os.system(f"xrdcp {inFile} .")
        fileName = inFile.split("/")[-1]
        skimmer = fileSkimmer(fileName,samps,"cuts.py",mode=samps['type'])
        skimmer.skim()
        os.system(f"xrdcp -f {skimmer.outFileName} root://cmseos.fnal.gov/{outDir}/")
        os.system(f"rm {skimmer.outFileName}")
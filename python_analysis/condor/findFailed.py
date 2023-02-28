import os
import shutil
import sys
import json
from XRootD import client
import numpy as np
import subprocess

if len(sys.argv) < 2:
    print("Wrong inputs!")
    print("Usage: python findFailed.py job_base_directory mode")

path = sys.argv[1]
mode = sys.argv[2]
if path[-1] == "/":
    path = path[:-1]
jobname_base = path.split("/")[-1]
dirs = [f for f in os.listdir(path) if os.path.isdir(path+"/"+f)]

xrdClient = client.FileSystem("root://cmseos.fnal.gov")
xrdPath = f"/store/group/lpcmetx/iDMe//analysis_output/{mode}/"+jobname_base+"/"
existing = [f.name for f in xrdClient.dirlist(xrdPath)[1] if ".coffea" in f.name]

for d in dirs:
    if f"output_{d}.coffea" not in existing:
        print(d)
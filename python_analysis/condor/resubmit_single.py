import os
import shutil
import sys
import json
from XRootD import client
import numpy as np
import subprocess

if len(sys.argv) < 2:
    print("Wrong inputs!")
    print("Usage: python resubmit_single.py submission_directory n_cores")

path = sys.argv[1]
n_cores = int(sys.argv[2])
if path[-1] == "/":
    path = path[:-1]
jobname = path.split("/")[-1]
jobname_base = path.split("/")[-2]

# re-make directories
tgz = "submissions/"+jobname_base+"/"+jobname+".tar.gz"
dirname = "submissions/"+jobname_base+"/"+jobname+"/"
if os.path.exists(tgz):
    os.remove(tgz)
if os.path.isdir(dirname+"Logs"):
    shutil.rmtree(dirname+"Logs")

# copy new versions of the analysis code
os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/analysisTools.py {0}".format(dirname))
os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/analysisSubroutines.py {0}".format(dirname))
os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/mySchema.py {0}".format(dirname))
os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/configs/histo_configs/histobins.py {0}".format(dirname))

subprocess.run("tar czf {0} -C submissions/{1}/{2}/ .".format(tgz,jobname_base,jobname),shell=True)
os.mkdir(dirname+"Logs/")

with open(dirname+"samples.json","r") as f:
    js = json.load(f)
    mode = js[0]['type']

submit_cmd = f"condor_submit condor_analysis.jdl -append \"Arguments = {jobname_base} {jobname} {mode} {n_cores}\" -append \"transfer_input_files = {tgz}\" -append \"Output = submissions/{jobname_base}/{jobname}/Logs/stdout.out\" -append \"Error = submissions/{jobname_base}/{jobname}/Logs/stderr.err\" -append \"Log = submissions/{jobname_base}/{jobname}/Logs/log.log\" -append \"request_cpus = {n_cores}\""

subprocess.run(submit_cmd,shell=True)
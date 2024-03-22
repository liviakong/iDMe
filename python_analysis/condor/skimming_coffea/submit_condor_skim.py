import os
import shutil
import sys
import json
from XRootD import client
import numpy as np
import subprocess

if len(sys.argv) < 5:
    print("Wrong inputs!")
    print("Usage: python submit_condor.py sampleConfig.json cutConfig.py n_file_per n_cores")

samples = sys.argv[1]
cuts = sys.argv[2]
n_file_per = int(sys.argv[3])
n_cores = int(sys.argv[4])

sampFile = samples.split("/")[-1]
cutFile = cuts.split("/")[-1]
jobname_base = sampFile.split(".")[0] + "_" + cutFile.split(".")[0]

xrdClient = client.FileSystem("root://cmseos.fnal.gov")

with open(samples) as f:
    sampData = json.load(f)

mode = sampData[0]["type"]
n_samp = len(sampData)
nfiles_by_samp = [samp['nFiles'] for samp in sampData]

for i in range(n_samp):
    samp = sampData[i]
    name = samp['name']
    outDir = "/store/group/lpcmetx/iDMe/skimmed_ntuples/{0}/{1}/output_{2}/".format(mode,jobname_base,name.replace(".","p"))
    print(outDir)
    xrdClient.mkdir(outDir,flags=client.flags.MkDirFlags.MAKEPATH)
    loc = samp['location']
    if type(loc) != list:
        fileList = ["root://cmsxrootd.fnal.gov/"+loc+f.name for f in xrdClient.dirlist(loc)[1] if '.root' in f.name]
    else:
        fileList = []
        for l in loc:
            fileList.extend(["root://cmsxrootd.fnal.gov/"+l+f.name for f in xrdClient.dirlist(l)[1] if '.root' in f.name])
    fileSets = [list(a) for a in np.array_split(fileList,1+len(fileList)//n_file_per)]
    job_idx = 0
    for i,fileSet in enumerate(fileSets):
        jobname = "ntuples_{0}_{1}".format(name,job_idx)
        dirname = "submissions_skim/"+jobname_base+"/"+jobname+"/"
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
        
        subsample = samp.copy()
        subsample["fileset"] = fileSet
        with open(dirname+"samples.json","w") as f:
            json.dump(subsample,f,indent=4)

        os.system("cp condor_skim.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/analysisTools.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/analysisSubroutines.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools/mySchema.py {0}".format(dirname))
        os.system("cp {0} {1}/cuts.py".format(cuts,dirname))

        tgz = "submissions_skim/"+jobname_base+"/"+jobname+".tar.gz"
        if os.path.exists(tgz):
            os.remove(tgz)
        subprocess.run("tar czf {0} -C submissions_skim/{1}/{2}/ .".format(tgz,jobname_base,jobname),shell=True)
        os.mkdir(f"submissions_skim/{jobname_base}/{jobname}/Logs/")

        print("done tarring")

        submit_cmd = f"condor_submit condor_skim.jdl -append \"Arguments = {jobname_base} {jobname} {mode} {n_cores} {outDir}\" -append \"transfer_input_files = {tgz}\" -append \"Output = submissions_skim/{jobname_base}/{jobname}/Logs/stdout.out\" -append \"Error = submissions_skim/{jobname_base}/{jobname}/Logs/stderr.err\" -append \"Log = submissions_skim/{jobname_base}/{jobname}/Logs/log.log\" -append \"request_cpus = {n_cores}\""
        
        subprocess.run(submit_cmd,shell=True)
        
        job_idx += 1
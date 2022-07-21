import os
import shutil
import sys
import json
from XRootD import client

if len(sys.argv) < 5:
    print("Wrong inputs!")
    print("Usage: python submit_condor.py sampleConfig.json histoConfig.json cutConfig.py n_files_per")

samples = sys.argv[1]
histos = sys.argv[2]
cuts = sys.argv[3]
n_file_per = int(sys.argv[4])

sampFile = samples.split("/")[-1]
histFile = histos.split("/")[-1]
jobname_base = sampFile.split(".")[0] + "_" + histFile.split(".")[0]

xrdClient = client.FileSystem("root://cmseos.fnal.gov")

with open(samples) as f:
    sampData = json.load(f)

mode = sampData[0]["type"]
n_samp = len(sampData)

outDir = "/store/group/lpcmetx/iDMe/analysis_output/{0}/{1}/".format(mode,jobname_base)
xrdClient.mkdir(outDir,flags=client.flags.MkDirFlags.MAKEPATH)

for i in range(n_samp):
    samp = sampData[i]
    name = samp["name"]
    nFiles = samp["nFiles"]
    fileList = ["root://cmsxrootd.fnal.gov/"+samp["location"]+f.name for f in xrdClient.dirlist(samp["location"])[1] if '.root' in f.name]
    job_idx = 0
    for i in range(0,len(fileList),n_file_per):
        fileSet = fileList[i:i+n_file_per]

        jobname = "{0}_{1}".format(name,job_idx)
        dirname = "submissions/"+jobname_base+"/"+jobname+"/"
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
        
        subsample = samp.copy()
        subsample["fileset"] = fileSet
        with open(dirname+"samples.json","w") as f:
            json.dump([subsample],f,indent=4)

        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/analysisTools/analysisTools.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/analysisTools/analysisSubroutines.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/analysisTools/mySchema.py {0}".format(dirname))
        os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/configs/histo_configs/histobins.py {0}".format(dirname))
        os.system("cp {0} {1}/histos.py".format(histos,dirname))
        os.system("cp {0} {1}/cuts.py".format(cuts,dirname))

        tgz = "submissions/"+jobname_base+"/"+jobname+".tar.gz"
        if os.path.exists(tgz):
            os.remove(tgz)
        os.system("tar czf {0} -C submissions/{1}/{2}/ .".format(tgz,jobname_base,jobname))

        submit_cmd = "condor_submit condor_analysis.jdl -append \"Arguments = {0} {1} {2}\" -append \"transfer_input_files = {3}\"".format(jobname_base,jobname,mode,tgz)
        os.system(submit_cmd)

        job_idx += 1
import os
import shutil
import sys
import json

if len(sys.argv) < 3:
    print("Wrong inputs!")
    print("Usage: python submit_condor.py sampleConfig.json histoConfig.json")

samples = sys.argv[1]
histos = sys.argv[2]

sampFile = samples.split("/")[-1]
histFile = histos.split("/")[-1]
jobname = sampFile.split(".")[0] + "_" + histFile.split(".")[0]

with open(samples) as f:
    sampData = json.load(f)

n_samp = len(sampData)
n_per = 5

for i in range(0,n_samp,n_per):
    idx = i // n_per
    dirname = "submissions/{0}_{1}/".format(jobname,idx)
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
    os.mkdir(dirname)
    subsample = sampData[i:i+n_per]
    with open(dirname+"samples.json","w") as f:
        json.dump(subsample,f,indent=4)
    
    os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/analysisTools/analysisTools.py {0}".format(dirname))
    os.system("cp /uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMeAnalysis/python_analysis/analysisTools/analysisSubroutines.py {0}".format(dirname))
    os.system("cp {0} {1}/histos.json".format(histos,dirname))

    if os.path.exists("submissions/{0}_{1}.tar.gz".format(jobname,idx)):
        os.remove("submissions/{0}_{1}.tar.gz".format(jobname,idx))
    os.system("tar czf submissions/{0}_{1}.tar.gz -C submissions/{0}_{1}/ .".format(jobname,idx))

    submit_cmd = "condor_submit condor_analysis.jdl -append \"Arguments = {0}_{1}\" -append \"transfer_input_files = submissions/{0}_{1}.tar.gz\"".format(jobname,idx)
    os.system(submit_cmd)
#!/usr/bin/env python
import coffea.util as util
import os
import json
import sys

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Bad input!")
        print("Usage: ./condor_analyze jobname_base jobname mode n_cores useGenAnalyzer")

    jobname_base = sys.argv[1]
    jobname = sys.argv[2]
    mode = sys.argv[3]
    n_cores = int(sys.argv[4])
    genAnalyzer = bool(sys.argv[5])

    os.system("tar xzf {0}.tar.gz".format(jobname))
    from analysisTools import Analyzer
    files = "samples.json"
    histos = "histos.py"
    cuts = "cuts.py"
    az = Analyzer(files,histos,cuts)
    if genAnalyzer:
        out = az.process(execr="futures",workers=n_cores,gen=True)
    else:
        out = az.process(execr="futures",workers=n_cores,gen=False)
    
    outName = "output_{0}.coffea".format(jobname)
    util.save(out,outName)

    outDir = "/store/group/lpcmetx/iDMe/analysis_output/{0}/{1}/".format(mode,jobname_base)

    copy_cmd = "xrdcp -f {0} root://cmseos.fnal.gov/{1}/{2}".format(outName,outDir,outName)
    os.system(copy_cmd)
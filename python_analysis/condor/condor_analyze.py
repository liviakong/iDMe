#!/usr/bin/env python
import coffea.util as util
import os
import json
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Bad input!")
        print("Usage: ./condor_analyze jobname")

    jobname = sys.argv[1]
    os.system("tar xzf {0}.tar.gz".format(jobname))
    from analysisTools import Analyzer
    files = "samples.json"
    histos = "histos.json"
    cuts = "cuts.py"
    az = Analyzer(files,histos,cuts)
    out = az.process(execr="futures")
    
    outName = "output_{0}.coffea".format(jobname)
    util.save(out,outName)

    copy_cmd = "xrdcp -f {0} root://cmseos.fnal.gov//store/group/lpcmetx/iDMe/analysis_output/signal/{0}".format(outName)
    os.system(copy_cmd)

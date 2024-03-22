import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("ERROR! Need 5 arguments")
        print("Usage: python3 submitCondor.py path_to_gridpack lifetime year nEvents nThreads [nev_per_job]")
        sys.exit()
    
    gridpack = sys.argv[1]
    lifetime = float(sys.argv[2])
    year = sys.argv[3]
    nevents = int(sys.argv[4])
    nThreads = int(sys.argv[5])

    hbar_c = 0.1973269804e-12 # in GeV * mm
    lifetime_width = hbar_c/lifetime

    procName = gridpack.split("/")[-1].split(".")[0]
    gridpackFileName = gridpack.split("/")[-1]
    baseDir = "submissions/"+procName
    here = os.getcwd()

    if not os.path.isdir(baseDir):
        os.mkdir(baseDir)
        os.mkdir(baseDir+"/Logs")
    if os.path.isdir(f"{baseDir}/submit"):
        os.system(f"rm -rf {baseDir}/submit")
    os.mkdir(f"{baseDir}/submit")
    os.system("cp genFragments/iDMe_pythiaGenFragment.py {0}".format(baseDir+"/submit"))
    os.system("cp genFromGridpack_*.sh {0}".format(baseDir+"/submit"))
    os.system("cp template_DIGIPremix_cfg_UL*.py {0}".format(baseDir+"/submit"))
    os.system("cp runGeneration.sh {0}".format(baseDir+"/submit"))
    os.system("cp {0} {1}".format(gridpack,baseDir+"/submit"))

    if os.path.exists(f"{baseDir}/submit.tar.gz"):
        os.system(f"rm {baseDir}/submit.tar.gz")
    os.chdir(baseDir+"/submit")
    os.system("tar -czf ../submit.tar.gz *")
    os.chdir(here)

    # make directory for log files
    logDir = baseDir+"/Logs/ctau-{0}_nev-{1}_year-{2}".format(lifetime,nevents,year)
    if os.path.isdir(logDir):
        os.system("rm -rf {0}".format(logDir))
    os.mkdir(logDir)

    if len(sys.argv) == 7:
        nev_per_job = int(sys.argv[6])
    else:
        nev_per_job = 10000
    print(f"{nev_per_job} events per job")
    nJobs = nevents // nev_per_job
    if nevents < nev_per_job:
        nJobs = 1
        nev_per_job = nevents
    
    condor_cmd = "condor_submit condorTemplate.jdl" 
    condor_cmd += " -append \"Arguments = {0} {1:.9e} {2} {3} {4} {5}\"".format(gridpackFileName,lifetime_width,lifetime,year,nev_per_job,nThreads)
    condor_cmd += " -append \"transfer_input_files = {0}/submit.tar.gz\"".format(baseDir)
    condor_cmd += " -append \"request_cpus = {0}\"".format(nThreads)
    condor_cmd += " -append \"output = {0}/\$(Cluster)_\$(Process).out\"".format(logDir)
    condor_cmd += " -append \"error = {0}/\\$(Cluster)_\$(Process).err\"".format(logDir)
    condor_cmd += " -append \"log = {0}/\\$(Cluster)_\\$(Process).log\"".format(logDir)
    condor_cmd += " -append \"Queue {0}\"".format(nJobs)

    print(condor_cmd)
    os.system(condor_cmd)


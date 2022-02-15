#!/usr/bin/env python

import os, sys

'''Usage: ./submit.py <LHE/gridpack filename> year [njobs]
'''

def buildExec(workpath, year):
    '''Given the workpath, write a exec.sh in it, to be used by condor'''

    execF = '''#!/bin/bash

export HOME=${PWD}

tar xvaf submit.tgz
cd submit
sh %s.sh $1
cd ${HOME}
rm -r submit/

exit 0'''
    with open(workpath + '/exec.sh', 'w') as f:
        f.write(execF % ('reprocessSignalNoTimingCleaning%s' % year))

def buildSubmit(workpath, uid, year):
    '''A series of actions to prepare submit dir'''

    stageOutPiece = '''
remoteDIR="/store/group/lpcmetx/iDMe/Samples/"
for f in `ls *AOD*.root`; do
    cmd="xrdcp -vf $f root://cmseos.fnal.gov/$remoteDIR/$f"
    echo $cmd && eval $cmd
done
'''
    os.mkdir(workpath+'/submit')
    os.system('cp reprocessSignalNoTimingCleaning%s.sh %s/submit/' % (year,workpath))
    #os.system('cp /uscms/home/sbrightt/nobackup/iDMe/CMSSW_patch_removeTimingCleaning/cmssw_10_6_26_noTimingClean.tar.gz %s/submit/' % workpath)
    #os.system('cp -r /uscms/home/sbrightt/nobackup/iDMe/CMSSW_patch_removeTimingCleaning/CMSSW_10_6_26 %s/submit/' % workpath)
    with open('%s/submit/reprocessSignalNoTimingCleaning%s.sh' % (workpath,year), 'a') as f:
        f.write(stageOutPiece)
    
    print "Tarring up submit..."
    os.chdir(workpath)
    os.system('tar -chzf submit.tgz submit')
    os.chdir('..')

def buildCondor(workpath, logpath, uid, year):
    '''build the condor file, return the abs path'''

    condorF = '''universe = vanilla
executable = {0}/exec.sh
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {0}/submit.tgz
transfer_output_files = ""
input = /dev/null
output = {1}/$(Cluster)_$(Process).out
error = {1}/$(Cluster)_$(Process).err
log = {1}/$(Cluster)_$(Process).log
rank = Mips
request_memory = 4000
+AcctGroup = "analysis"
+ProjectName = "DarkMatterSimulation"
Queue 1'''.format(workpath, logpath)
    condorFN = 'condor_config.jdl'

    with open(workpath + '/' + condorFN, 'w') as jdlfile:
        jdlfile.write(condorF)

    return os.path.join(workpath, condorFN)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "ERROR! Need at least 2 arguments!"
        print "Usage: ./runReprocessNoTiming.py fileList year"
        sys.exit()

    year = sys.argv[2]
    flist = sys.argv[1]
    listName = flist.split("/")[-1].split(".")[0]
    files = []
    with open(flist) as f:
        files = f.read().splitlines()

    Logpath = os.getcwd() + '/Logs'
    Submissionpath = os.getcwd() + '/submissions'
    if not os.path.isdir(Logpath): os.mkdir(Logpath)
    if not os.path.isdir(Submissionpath): os.mkdir(Submissionpath)
    Workpath = Submissionpath + '/submit_reprocessNoTiming_'+listName
    if os.path.isdir(Workpath): os.system('rm -rf %s' % Workpath)
    os.mkdir(Workpath)
    Uid = os.getuid()

    buildSubmit(workpath=Workpath, uid=Uid, year=year)
    buildExec(workpath=Workpath,year=year)
    theCondor = buildCondor(workpath=Workpath,
                logpath=Logpath, uid=Uid, year=year)
    for f in files:
        submit_string = theCondor + " -append \"Arguments = {0}\"".format(f)
        #print(submit_string)
        os.system('condor_submit %s' % submit_string)

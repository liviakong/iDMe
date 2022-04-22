from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
import CRABClient.UserUtilities
from optparse import OptionParser
from multiprocessing import Process
import os
import json 

def parseArguments():
    parser = OptionParser()
    
    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-y', '--year',
                      dest = 'year',
                      default = '2018',
                      help = "Which year to process ('2018'(default)/'2017'/'2016')",
                      metavar = 'YEAR')

    parser.add_option("-f","--inFile",
                      dest = 'inFile',
                      default = '',
                      help = 'JSON file with input data or MC lists',
                      metavar = 'INFILE')
    
    (options, arguments) = parser.parse_args()
    
    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if not options.year:
        parser.error("-y year option not provided")
    if not options.inFile:
        parser.error("-f input file option not provided")

    return options

def main():
    options = parseArguments()
    
    if 'CMSSW_BASE' not in os.environ.keys():
        print "Run cmsenv first!"
        return
    base_dir = os.environ['CMSSW_BASE']

    year = options.year

    config = CRABClient.UserUtilities.config()

    # Basic settings common to all runs 
    config.General.workArea = base_dir+'/src/iDMeAnalysis/iDMeSkimmer/crab/'
    config.General.transferOutputs = True
    config.General.transferLogs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = base_dir+'/src/iDMeAnalysis/iDMeSkimmer/scripts/run_ntuplizer_cfg.py'
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.numCores = 1
    #config.Data.splitting = 'Automatic'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 30
    config.Data.publication = False
    config.Data.ignoreLocality = True
    config.Site.whitelist = ['T2_US_*', 'T2_DE_*', 'T2_EE_*', 'T2_ES_*', 'T2_GR_*', 'T2_HU_*', 'T2_IT_*', 'T2_RU_*', 'T2_UK_*']
    config.Site.storageSite = 'T3_US_FNALLPC'

    inFile = options.inFile
    samples = {}
    with open(inFile) as f:
        samples = json.load(f)
    for samp in samples.keys():
        output_base = '/store/group/lpcmetx/iDMe/Samples/Ntuples/background/{0}/'.format(year)
        config.Data.outLFNDirBase = output_base
        for subsample, dataset in samples[samp].items():
            config.Data.inputDataset = dataset
            config.General.requestName = 'iDMe_'+subsample
            config.JobType.pyCfgParams = ['numThreads=1','outfile={0}.root'.format(subsample)]
            print 'Submitting for input dataset {0}'.format(subsample)
            #crabCommand(options.crabCmd, config = config)
            kwargs = {'config':config}
            p = Process(target=crabCommand,args=(options.crabCmd,),kwargs=kwargs)
            p.start()
            p.join()

if __name__ == '__main__':
    main()

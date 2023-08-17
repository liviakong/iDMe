from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_ntuplizer_test0'
config.General.workArea = '/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/iDMeSkimmer/crab/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/iDMeSkimmer/scripts/run_ntuplizer_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['outfile=test_QCD.root','numThreads=1']
config.JobType.numCores = 1

config.section_("Data")
config.Data.splitting = 'Automatic'
config.Data.outLFNDirBase = '/store/group/lpcmetx/iDMe/Samples/Ntuples/background/QCD/'
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.inputDataset = '/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'

config.section_("Site")
config.Site.whitelist = ['T2_US_*', 'T2_DE_*', 'T2_EE_*', 'T2_ES_*', 'T2_GR_*', 'T2_HU_*', 'T2_IT_*', 'T2_RU_*', 'T2_UK_*']
config.Site.storageSite = 'T3_US_FNALLPC'

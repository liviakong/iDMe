# Trying to frankenstein together the miniAOD producer and the ntuplizer

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
import json

# argument parsing
options = VarParsing.VarParsing('analysis')
options.register('data',
        False,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.bool,
        "Run on data (1) or MC (0)"
        )
options.register('signal',
        True,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.bool,
        "Run on signal (1) or not (0")
options.register('year',
        2018,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
        "Data/MC year")
options.register('numThreads',
        8,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
        "Number of threads (for CRAB vs non-CRAB execution)")
options.register("nEvents",
	-1,
	VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
	"Number of events to process (defaults to all)")
options.register('flist',
        "",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File list to ntuplize")
options.register('outfile',
        "test_output.root",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "Output file name")

options.parseArguments()

# file list
if options.flist != "":
    if ".txt" in options.flist:
        # list of files
        print("reading input file list: "+options.flist)
        options.inputFiles = FileUtils.loadListFromFile(options.flist)
    else:
        # we have passed a file name directly
        options.inputFiles = options.flist

# globaltag
globaltag = ''
if options.year == 2016:
    globaltag = '106X_mcRun2_asymptotic_v13'
    era = Run2_2016
elif options.year == 2017:
    globaltag = '106X_mc2017_realistic_v6'
    era = Run2_2017
elif options.year == 2018:
    globaltag = '106X_upgrade2018_realistic_v16_L1v1'
    era = Run2_2018
else:
    print("Invalid year given for run 2 : {0}".format(options.year))
    exit

# Define the process
process = cms.Process('PAT',era,run2_miniAOD_UL)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Input source
if options.inputFiles is not None:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
        skipBadFiles = cms.untracked.bool(True)
    )
else:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring(),
        skipBadFiles = cms.untracked.bool(True)
    )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:{0}'.format(options.nEvents)),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Setting the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

# max events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
)

# Output definition

"""process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:output.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)"""
###########################################
###### Path and EndPath definitions #######
###########################################

# Various filters (MiniAOD producer)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_hfNoisyHitsFilter = cms.Path(process.hfNoisyHitsFilter)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
# MiniAOD end of job steps
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput) ### we don't need the MiniAOD file to actually be written!

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,
                                process.Flag_HBHENoiseIsoFilter,
                                process.Flag_CSCTightHaloFilter,
                                process.Flag_CSCTightHaloTrkMuUnvetoFilter,
                                process.Flag_CSCTightHalo2015Filter,
                                process.Flag_globalTightHalo2016Filter,
                                process.Flag_globalSuperTightHalo2016Filter,
                                process.Flag_HcalStripHaloFilter,
                                process.Flag_hcalLaserEventFilter,
                                process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                process.Flag_EcalDeadCellBoundaryEnergyFilter,
                                process.Flag_ecalBadCalibFilter,
                                process.Flag_goodVertices,
                                process.Flag_eeBadScFilter,
                                process.Flag_ecalLaserCorrFilter,
                                process.Flag_trkPOGFilters,
                                process.Flag_chargedHadronTrackResolutionFilter,
                                process.Flag_muonBadTrackFilter,
                                process.Flag_BadChargedCandidateFilter,
                                process.Flag_BadPFMuonFilter,
                                process.Flag_BadPFMuonDzFilter,
                                process.Flag_hfNoisyHitsFilter,
                                process.Flag_BadChargedCandidateSummer16Filter,
                                process.Flag_BadPFMuonSummer16Filter,
                                process.Flag_trkPOG_manystripclus53X,
                                process.Flag_trkPOG_toomanystripclus53X,
                                process.Flag_trkPOG_logErrorTooManyClusters,
                                process.Flag_METFilters,
                                process.endjob_step)
                                #process.MINIAODSIMoutput_step)

# Associate the pat tasks
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.numThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# Adding ntuplizer process as a subprocess
ntupleProcess = cms.Process("USER",era,run2_miniAOD_UL)
process.addSubProcess(cms.SubProcess(ntupleProcess))

ntupleProcess.load("FWCore.MessageService.MessageLogger_cfi")
ntupleProcess.load('Configuration.StandardSequences.Services_cff')
ntupleProcess.load("Configuration.EventContent.EventContent_cff")
ntupleProcess.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
ntupleProcess.load("Configuration.StandardSequences.GeometryRecoDB_cff")
ntupleProcess.load('Configuration.StandardSequences.MagneticField_38T_cff')
ntupleProcess.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

ntupleProcess.GlobalTag.globaltag = globaltag

ntupleProcess.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.numThreads)
    )
ntupleProcess.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
    )
ntupleProcess.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outfile),
    closeFileFast = cms.untracked.bool(True)
    )

#######################
##### MET Filters #####
#######################
metFilters = []
if options.year == 2016:
    metFilters = [
        "Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadPFMuonDzFilter",
        "Flag_eeBadScFilter",
        "Flag_hfNoisyHitsFilter"
    ]
elif options.year == 2017 or options.year == 2018:
    metFilters = [
        "Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadPFMuonDzFilter",
        "Flag_hfNoisyHitsFilter",
        "Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter"
    ]

#######################
###### Triggers #######
#######################
triggerPaths16 = [
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250",
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele25_WPTight_Gsf",
    "HLT_Ele27_WPTight_Gsf",
    "HLT_PFMET90_PFMHT90_IDTight",
    "HLT_PFMET100_PFMHT100_IDTight",
    "HLT_PFMET110_PFMHT110_IDTight",
    "HLT_PFMET120_PFMHT120_IDTight"
]

triggerPaths17 = [
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele27_WPTight_Gsf",
    "HLT_Ele32_WPTight_Gsf",
    "HLT_Ele35_WPTight_Gsf",
    "HLT_Ele38_WPTight_Gsf",
    "HLT_Ele40_WPTight_Gsf",
    "HLT_PFMET110_PFMHT110_IDTight",
    "HLT_PFMET120_PFMHT120_IDTight",
    "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight",
    "HLT_PFMETTypeOne110_PFMHT110_IDTight",
    "HLT_PFMETTypeOne120_PFMHT120_IDTight",
    "HLT_PFMETTypeOne130_PFMHT130_IDTight",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight"
]

triggerPaths18 = [
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele27_WPTight_Gsf",
    "HLT_Ele28_WPTight_Gsf",
    "HLT_Ele30_WPTight_Gsf",
    "HLT_Ele32_WPTight_Gsf",
    "HLT_Ele35_WPTight_Gsf",
    "HLT_Ele38_WPTight_Gsf",
    "HLT_Ele40_WPTight_Gsf",
    "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
    "HLT_PFMET110_PFMHT110_IDTight",
    "HLT_PFMET120_PFMHT120_IDTight",
    "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight",
    "HLT_PFMETTypeOne110_PFMHT110_IDTight",
    "HLT_PFMETTypeOne120_PFMHT120_IDTight",
    "HLT_PFMETTypeOne130_PFMHT130_IDTight",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight"
]

# Electron effective area input file for PU-corrected PF isolation calculations
effAreaInputPath = "RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"

##############################
###### Main iDM analyzer #####
##############################
## (note: cfi with default options is generated by fillDescriptions when compiling -- not found in python/)
## (the only options specified below after cloning are the non-default ones)
import iDMe
print(dir(iDMe))
from iDMe.AODSkimmer.ElectronSkimmer_cfi import ElectronSkimmer
# setting up skimmer & ntupleProcess path
ntupleProcess.ntuples = ElectronSkimmer.clone(
    isData = cms.bool(options.data),
    isSignal = cms.bool(options.signal),
    year=options.year,
    metFilters = cms.vstring(metFilters),
    triggerPaths16 = cms.vstring(triggerPaths16),
    triggerPaths17 = cms.vstring(triggerPaths17),
    triggerPaths18 = cms.vstring(triggerPaths18),
    allTriggerPaths = cms.vstring(list(set(triggerPaths16+triggerPaths17+triggerPaths18))),
    effAreasConfigFile = cms.FileInPath(effAreaInputPath)
)
ntupleProcess.totalPath = cms.Path(cms.Sequence(ntupleProcess.ntuples))
ntupleProcess.schedule = cms.Schedule(ntupleProcess.totalPath)
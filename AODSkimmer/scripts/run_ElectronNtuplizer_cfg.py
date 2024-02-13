import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
import json

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

# reading input file list if given
if options.flist != "":
    if ".txt" in options.flist:
        # list of files
        print("reading input file list: "+options.flist)
        options.inputFiles = FileUtils.loadListFromFile(options.flist)
    else:
        # we have passed a file name directly
        options.inputFiles = options.flist 
        print("running on file:")
        print(options.inputFiles)

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

process = cms.Process("USER",era,run2_miniAOD_UL)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.EventContent.EventContent_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = globaltag

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.numThreads)
    )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
    )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    skipBadFiles = cms.untracked.bool(True)
    )
process.TFileService = cms.Service("TFileService",
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
from iDMe.AODSkimmer.ElectronSkimmer_cfi import ElectronSkimmer
# setting up skimmer & process path
process.ntuples = ElectronSkimmer.clone(
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
process.commonSequence = cms.Sequence(
    process.ntuples
)

if options.data:
    process.totalPath = cms.Path(
        process.commonSequence
    )
else:
    process.totalPath = cms.Path(
        process.commonSequence
    )

process.schedule = cms.Schedule(process.totalPath)

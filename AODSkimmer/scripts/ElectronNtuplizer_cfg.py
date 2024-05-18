import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.Eras.Era_Run2_2016_HIPM_cff import Run2_2016_HIPM
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
import json
import sys

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
        "2018",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
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
if ".txt" in options.flist:
    # list of files
    print("reading input file list: "+options.flist)
    options.inputFiles = FileUtils.loadListFromFile(options.flist)
else:
    # we have passed a file name directly
    options.inputFiles = options.flist

# globaltag
globaltag = ''
if options.year == '2016APV':
    globaltag = '106X_dataRun2_v37' if options.data else '106X_mcRun2_asymptotic_preVFP_v11'
    era = Run2_2016_HIPM
    recoEgammaTools_era = '2016preVFP-UL'
elif options.year == '2016':
    globaltag = '106X_dataRun2_v37' if options.data else '106X_mcRun2_asymptotic_v17'
    era = Run2_2016
    recoEgammaTools_era = '2016postVFP-UL'
elif options.year == '2017':
    globaltag = '106X_dataRun2_v37' if options.data else '106X_mc2017_realistic_v10'
    era = Run2_2017
    recoEgammaTools_era = '2017-UL'
elif options.year == '2018':
    globaltag = '106X_dataRun2_v37' if options.data else '106X_upgrade2018_realistic_v16_L1v1'
    era = Run2_2018
    recoEgammaTools_era = '2018-UL'
else:
    print("Invalid year given for run 2 : {0}".format(options.year))
    exit

#######################
##### MET Filters #####
#######################
metFilters = []
if options.year == '2016' or options.year == '2016APV':
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
elif options.year == '2017' or options.year == '2018':
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
# record all trigger paths that might be useful acrcoss all years - some will not always be available,
# but what's available will get written out to the ntuples
triggerPaths = [
    # MET triggers
    "HLT_PFMET90_PFMHT90_IDTight",
    "HLT_PFMET100_PFMHT100_IDTight",
    "HLT_PFMET110_PFMHT110_IDTight",
    "HLT_PFMET120_PFMHT120_IDTight",
    "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight",
    "HLT_PFMETTypeOne110_PFMHT110_IDTight",
    "HLT_PFMETTypeOne120_PFMHT120_IDTight",
    "HLT_PFMETTypeOne130_PFMHT130_IDTight",
    "HLT_PFMETTypeOne140_PFMHT140_IDTight",
    "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v9"
    # Jet triggers (for MET trigger eff)
    "HLT_PFJet15",
    "HLT_PFJet25",
    "HLT_PFJet40",
    "HLT_PFJet60",
    "HLT_PFJet80",
    "HLT_PFJet140",
    "HLT_PFJet200",
    "HLT_PFJet260",
    "HLT_PFJet320",
    "HLT_PFJet400",
    "HLT_PFJet450",
    "HLT_PFJet500",
    "HLT_PFJet550",
    "HLT_AK4PFJet30",
    "HLT_AK4PFJet50",
    "HLT_AK4PFJet80",
    "HLT_AK4PFJet100",
    "HLT_AK4PFJet120",
    # Electron triggers (for MET trigger eff)
    "HLT_DoubleEle25_CaloIdL_MW",
    "HLT_DoubleEle27_CaloIdL_MW",
    "HLT_DoubleEle33_CaloIdL_MW",
    "HLT_DoubleEle24_eta2p1_WPTight_Gsf",
    "HLT_Ele20_WPTight_Gsf",
    "HLT_Ele15_WPLoose_Gsf",
    "HLT_Ele17_WPLoose_Gsf",
    "HLT_Ele20_WPLoose_Gsf",
    "HLT_Ele27_WPTight_Gsf",
    "HLT_Ele28_WPTight_Gsf",
    "HLT_Ele30_WPTight_Gsf",
    "HLT_Ele32_WPTight_Gsf",
    "HLT_Ele35_WPTight_Gsf",
    "HLT_Ele38_WPTight_Gsf",
    "HLT_Ele40_WPTight_Gsf",
    "HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
    "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",
    "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
    "HLT_Ele17_CaloIdM_TrackIdM_PFJet30",
    "HLT_Ele23_CaloIdM_TrackIdM_PFJet30",
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"
]

# Electron effective area input file for PU-corrected PF isolation calculations
effAreaInputPath = "RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"

process = cms.Process("USER",era,run2_miniAOD_UL)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.EventContent.EventContent_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

process.options = cms.untracked.PSet()
process.options.numberOfThreads=cms.untracked.uint32(options.numThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

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

##############################
###### Main iDM analyzer #####
##############################
# import ntuplizer
from iDMe.AODSkimmer.ElectronSkimmer_cfi import ElectronSkimmer
process.ntuples = ElectronSkimmer.clone(
    isData = cms.bool(options.data),
    isSignal = cms.bool(options.signal),
    year = options.year,
    metFilters = cms.vstring(metFilters),
    triggerPaths = cms.vstring(triggerPaths),
    effAreasConfigFile = cms.FileInPath(effAreaInputPath)
)

# import EGamma postreco tools
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era=recoEgammaTools_era)

# load nanoAOD producer chain for low-pT electrons -- computes mini iso
process.load('PhysicsTools.NanoAOD.lowPtElectrons_cff')
process.lowPtNanoElectronSequence = cms.Sequence(process.modifiedLowPtElectrons+\
                                         process.updatedLowPtElectrons+\
                                         process.lowPtPATElectronID+\
                                         process.isoForLowPtEle+\
                                         process.updatedLowPtElectronsWithUserData)

# load nanoAOD producer for electrons - modify to just compute PFiso and miniIso
process.load('PhysicsTools.NanoAOD.electrons_cff')
process.isoForEleRelative = process.isoForEle.clone( # configure isolation to compute relative isolation by default
    relative = cms.bool(True)
)
process.slimmedElectronsWithUserDataMinimal = process.slimmedElectronsWithUserData.clone( # remove need to save non-isolation value maps
    src = cms.InputTag("slimmedElectrons"),
    userFloats = cms.PSet(
        miniIsoChg = cms.InputTag("isoForEleRelative:miniIsoChg"),
        miniIsoAll = cms.InputTag("isoForEleRelative:miniIsoAll"),
        PFIsoChg = cms.InputTag("isoForEleRelative:PFIsoChg"),
        PFIsoAll = cms.InputTag("isoForEleRelative:PFIsoAll"),
        PFIsoAll04 = cms.InputTag("isoForEleRelative:PFIsoAll04"),
    ),
    userIntFromBools = cms.PSet(),
    userInts = cms.PSet(),
    userCands = cms.PSet()
)
process.nanoElectronSequence = cms.Sequence(process.isoForEleRelative+\
                                            process.slimmedElectronsWithUserDataMinimal)

## Define paths and schedule
process.ntupleSequence = cms.Sequence(process.ntuples)
process.ntuplePath = cms.Path(process.ntupleSequence)

process.iDMEgammaPostRecoSequence = cms.Sequence(process.egammaPostRecoSeq)
process.iDMEgammaPostReco = cms.Path(process.iDMEgammaPostRecoSequence)

process.iDMNanoElectronSequence = cms.Sequence(process.lowPtNanoElectronSequence + process.nanoElectronSequence)
process.iDMNanoElectron = cms.Path(process.iDMNanoElectronSequence)

process.schedule = cms.Schedule(process.iDMEgammaPostReco,process.iDMNanoElectron,process.ntuplePath)

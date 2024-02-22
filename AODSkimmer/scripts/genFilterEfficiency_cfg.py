import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils

options = VarParsing.VarParsing('analysis')
options.register('flist',
        "",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File list to ntuplize")
options.register("nEvents",
	-1,
	VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
	"Number of events to process (defaults to all)")
options.register('outfile',
        "filter_eff.root",
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
        options.inputFiles = cms.untracked.vstring(options.flist)

process = cms.Process("GenFilterEfficiency")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outfile),
    closeFileFast = cms.untracked.bool(True)
)

#process.load("GeneratorInterface.Core")
#process.dummy = process.GenFilterEfficiencySaver.clone(
#    genFilterInfoTag = cms.InputTag("genFilterEfficiencyProducer")
#)
process.efficiency = cms.EDAnalyzer("GenFilterEfficiencySaver",
                               genFilterInfoTag = cms.InputTag("genFilterEfficiencyProducer")
)

process.p = cms.Path(process.efficiency)
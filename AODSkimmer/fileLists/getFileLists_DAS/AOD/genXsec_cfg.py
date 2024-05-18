import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.Utilities.FileUtils as FileUtils
options = VarParsing ('analysis')
options.parseArguments()
process = cms.Process('XSec')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

fileList = FileUtils.loadListFromFile(options.inputFiles[0])
fileList = ["root://cmsxrootd.fnal.gov/"+f for f in fileList]

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(fileList), 
    secondaryFileNames = secFiles)
process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")

process.ana = cms.Path(process.xsec)
process.schedule = cms.Schedule(process.ana)

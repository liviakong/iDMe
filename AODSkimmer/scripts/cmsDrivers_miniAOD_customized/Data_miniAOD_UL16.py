# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 -s PAT --runUnscheduled --data --era Run2_2016 --scenario pp --conditions 106X_dataRun2_v37 --eventcontent MINIAOD --datatier MINIAOD --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016,Configuration/DataProcessing/Utils.addMonitoring --filein test.root --python_filename=Data_miniAOD_UL16.py --no_exec
import FWCore.ParameterSet.Config as cms

def makeProcess(era,procModifiers,globalTag,options):
    process = cms.Process('PAT',era,*procModifiers)

    # import of standard configurations
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.load('Configuration.EventContent.EventContent_cff')
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
    process.load('Configuration.StandardSequences.PAT_cff')
    process.load('Configuration.StandardSequences.EndOfProcess_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(options.nEvents)
    )

    # Input source
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
        secondaryFileNames = cms.untracked.vstring(),
        skipBadFiles = cms.untracked.bool(True)
    )

    # custom options
    process.options = cms.untracked.PSet()
    process.options.numberOfThreads=cms.untracked.uint32(options.numThreads)
    process.options.numberOfStreams=cms.untracked.uint32(0)
    process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

    # TFileService for ntuple output
    process.TFileService = cms.Service("TFileService",
                                        fileName = cms.string(options.outfile),
                                        closeFileFast = cms.untracked.bool(True)
                                        )

    # Production Info
    process.configurationMetadata = cms.untracked.PSet(
        annotation = cms.untracked.string('step1 nevts:{0}'.format(options.nEvents)),
        name = cms.untracked.string('Applications'),
        version = cms.untracked.string('$Revision: 1.19 $')
    )

    # Other statements
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

    # Path and EndPath definitions
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
    process.endjob_step = cms.EndPath(process.endOfProcess)

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
    
    process.schedule.associate(process.patTask)
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)

    # customisation of the process.

    # Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
    from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2016 

    #call to customisation function customisePostEra_Run2_2016 imported from Configuration.DataProcessing.RecoTLR
    process = customisePostEra_Run2_2016(process)

    # Automatic addition of the customisation function from Configuration.DataProcessing.Utils
    from Configuration.DataProcessing.Utils import addMonitoring 

    #call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
    process = addMonitoring(process)

    # End of customisation functions
    #do not add changes to your config after this point (unless you know what you are doing)
    # don't want unscheduled execution - need to run my own paths!
    #from FWCore.ParameterSet.Utilities import convertToUnscheduled
    #process=convertToUnscheduled(process)

    # customisation of the process.

    # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

    #call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    process = miniAOD_customizeAllData(process)

    # End of customisation functions

    # Customisation from command line

    # Add early deletion of temporary data products to reduce peak memory need
    from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
    process = customiseEarlyDelete(process)
    # End adding early deletion

    return process
import FWCore.ParameterSet.Config as cms

# Loading MET filters recommended for UL (https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#UL_data)
# Setting up MET filters as in https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py

# Custom primary vertex filter that allows taggingMode = True
from iDMe.AODSkimmer.myPrimaryVertexFilter_cfi import *
# Beam halo filter
from RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi import *
# HBHE filters
from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
# ECAL dead cell primitive filters
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
# Bad PF muon filters
from RecoMET.METFilters.BadPFMuonFilter_cfi import *
from RecoMET.METFilters.BadPFMuonDzFilter_cfi import *
# HF noisy hits filter
from RecoMET.METFilters.hfNoisyHitsFilter_cfi import *
# eeBadSc filter
from RecoMET.METFilters.eeBadScFilter_cfi import *
# ECAL bad calibration filter
from RecoMET.METFilters.ecalBadCalibFilter_cfi import *

# Set tagging mode to True, so we save filter result rather than skipping event entirely
myPrimaryVertexFilter.taggingMode = cms.bool(True)
globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)
HBHENoiseFilter.taggingMode = cms.bool(True)
HBHENoiseIsoFilter.taggingMode = cms.bool(True)
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
BadPFMuonFilter.taggingMode = cms.bool(True)
BadPFMuonDzFilter.taggingMode = cms.bool(True)
hfNoisyHitsFilter.taggingMode = cms.bool(True)
eeBadScFilter.taggingMode = cms.bool(True)
ecalBadCalibFilter.taggingMode = cms.bool(True)

# creating filter sequences
myMetFilters = cms.Sequence(
    HBHENoiseFilterResultProducer*
    myPrimaryVertexFilter*
    globalSuperTightHalo2016Filter*
    HBHENoiseFilter*
    HBHENoiseIsoFilter*
    EcalDeadCellTriggerPrimitiveFilter*
    BadPFMuonFilter*
    BadPFMuonDzFilter*
    hfNoisyHitsFilter*
    eeBadScFilter*
    ecalBadCalibFilter
)

# exclude some filters from phase2_hgcal era (as in https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/RecoMET/METFilters/python/metFilters_cff.py)
from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
phase2_hgcal.toReplaceWith(myMetFilters, myMetFilters.copyAndExclude([
    HBHENoiseFilterResultProducer, HBHENoiseFilter, # No hcalnoise for hgcal
    eeBadScFilter                                   # No EE
]))
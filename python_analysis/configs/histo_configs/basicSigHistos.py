from histobins import *
from coffea.hist import Hist
import numpy as np
import awkward as ak

histograms = {
    "ele_kinematics" : Hist("Events",sample,cut,ele_type,ele_pt,ele_eta,ele_phi),
    "gen_ee_kinematics" : Hist("Events",sample,cut,mass,dR),
    "ele_trackHits" : Hist("Events",sample,cut,ele_type,ele_trkHits,ele_pixHits,ele_stripHits),
    "ele_trackQual" : Hist("Events",sample,cut,ele_type,ele_chi2,ele_trkIso),
    "gen_displacement" : Hist("Events",sample,cut,vxy,vz),
    "RRvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "RRvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
    "LLvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "LLvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
    "LRvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "LRvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
}

subroutines = []

def fillHistos(events,histos,samp,cut):
    # Electron information
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Default",
                pt=ak.flatten(events.Electron.pt),
                eta=ak.flatten(events.Electron.eta),
                phi=ak.flatten(events.Electron.phi))
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Low pT",
                pt=ak.flatten(events.LptElectron.pt),
                eta=ak.flatten(events.LptElectron.eta),
                phi=ak.flatten(events.LptElectron.phi))
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Candidate",
                pt=ak.flatten(events.EleCand.pt),
                eta=ak.flatten(events.EleCand.eta),
                phi=ak.flatten(events.EleCand.phi))
    
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Default",
                    numTrkHits=ak.flatten(events.Electron.numTrackerHits),
                    numPixHits=ak.flatten(events.Electron.numPixHits),
                    numStripHits=ak.flatten(events.Electron.numStripHits))
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Low pT",
                    numTrkHits=ak.flatten(events.LptElectron.numTrackerHits),
                    numPixHits=ak.flatten(events.LptElectron.numPixHits),
                    numStripHits=ak.flatten(events.LptElectron.numStripHits))
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Candidate",
                    numTrkHits=ak.flatten(events.EleCand.numTrackerHits),
                    numPixHits=ak.flatten(events.EleCand.numPixHits),
                    numStripHits=ak.flatten(events.EleCand.numStripHits))

    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Default",
                    chi2=ak.flatten(events.Electron.trkChi2),
                    trkIso=ak.flatten(events.Electron.trkIso))
    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Low pT",
                    chi2=ak.flatten(events.LptElectron.trkChi2),
                    trkIso=ak.flatten(events.LptElectron.trkIso))
    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Candidate",
                    chi2=ak.flatten(events.EleCand.trkChi2),
                    trkIso=ak.flatten(events.EleCand.trkIso))

    # Filling gen particle information (signal only)
    gen_eles = ak.concatenate([events.GenEle,events.GenPos])
    histos['ele_kinematics'].fill(sample=samp,cut=cut,ele_type="Generator",
                        pt=gen_eles.pt,eta=gen_eles.eta,phi=gen_eles.phi)
    histos['gen_displacement'].fill(sample=samp,cut=cut,
                        vxy=gen_eles.vxy,
                        vz=gen_eles.vz)
    histos["gen_ee_kinematics"].fill(sample=samp,cut=cut,
                mass=events.genEE.mass,
                dR=events.genEE.dr)

    # Reco vertex information
    histos["RRvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.RRvtx.vxy),
                        sigma_vxy=ak.flatten(events.RRvtx.sigmavxy))
    histos["RRvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.RRvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.RRvtx.vxy/events.RRvtx.sigmavxy))
    # low-low, i.e. two low-pT electrons matched to a vertex
    histos["LLvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.LLvtx.vxy),
                        sigma_vxy=ak.flatten(events.LLvtx.sigmavxy))
    histos["LLvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.LLvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.LLvtx.vxy/events.LLvtx.sigmavxy))
    # low-reg, i.e. one low-pT and one default electron matched to a vertex
    histos["LRvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.LRvtx.vxy),
                        sigma_vxy=ak.flatten(events.LRvtx.sigmavxy))
    histos["LRvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.LRvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.LRvtx.vxy/events.LRvtx.sigmavxy))


from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        # Selected electron 1 histos, 1D
        #"sel_e1_pt" : Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        #"sel_e1_trkIso" : Hist(samp,cut,ele_trkIso,storage=hist.storage.Weight()),
        #"sel_e1_trkRelIso" : Hist(samp,cut,ele_trkRelIso,storage=hist.storage.Weight()),
        #"sel_e1_PFRelIso" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e1_PFRelIso3" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e1_PFRelIso4" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e1_PFRelIso8" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e1_PFIso3" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e1_PFIso4" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e1_PFIso8" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e1_trkChi2" : Hist(samp,cut,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e1_trkProb" : Hist(samp,cut,ele_prob,storage=hist.storage.Weight()),
        #"sel_e1_dxy" : Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        #"sel_e1_dxySignif" : Hist(samp,cut,ele_dxySignif,storage=hist.storage.Weight()),
        #"sel_e1_angRes" : Hist(samp,cut,ele_angRes,storage=hist.storage.Weight()),
        #"sel_e1_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        #"sel_e1_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        # Selected electron 2 histos, 1D
        #"sel_e2_pt" : Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        #"sel_e2_trkIso" : Hist(samp,cut,ele_trkIso,storage=hist.storage.Weight()),
        #"sel_e2_trkRelIso" : Hist(samp,cut,ele_trkRelIso,storage=hist.storage.Weight()),
        #"sel_e2_PFRelIso" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e2_PFRelIso3" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e2_PFRelIso4" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e2_PFRelIso8" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_e2_PFIso3" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e2_PFIso4" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e2_PFIso8" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_e2_trkChi2" : Hist(samp,cut,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e2_trkProb" : Hist(samp,cut,ele_prob,storage=hist.storage.Weight()),
        #"sel_e2_dxy" : Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        #"sel_e2_dxySignif" : Hist(samp,cut,ele_dxySignif,storage=hist.storage.Weight()),
        #"sel_e2_angRes" : Hist(samp,cut,ele_angRes,storage=hist.storage.Weight()),
        #"sel_e2_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        #"sel_e2_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        # Selected electron 1 histos, 2D
        #"sel_e1_pt_vs_trkIso" : Hist(samp,cut,ele_pt,ele_trkIso,storage=hist.storage.Weight()),
        #"sel_e1_pt_vs_trkRelIso" : Hist(samp,cut,ele_pt,ele_trkRelIso,storage=hist.storage.Weight()),
        #"sel_e1_pt_vs_chi2" : Hist(samp,cut,ele_pt,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e1_chi2_vs_dxy" : Hist(samp,cut,ele_chi2,ele_dxy,storage=hist.storage.Weight()),
        #"sel_e1_numHits_vs_trkChi2" : Hist(samp,cut,ele_trkHits,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e1_numHits_vs_trkProb" : Hist(samp,cut,ele_trkHits,ele_prob,storage=hist.storage.Weight()),
        # Selected electron 2 histos, 2D
        #"sel_e2_pt_vs_trkIso" : Hist(samp,cut,ele_pt,ele_trkIso,storage=hist.storage.Weight()),
        #"sel_e2_pt_vs_trkRelIso" : Hist(samp,cut,ele_pt,ele_trkRelIso,storage=hist.storage.Weight()),
        #"sel_e2_pt_vs_chi2" : Hist(samp,cut,ele_pt,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e2_chi2_vs_dxy" : Hist(samp,cut,ele_chi2,ele_dxy,storage=hist.storage.Weight()),
        #"sel_e2_numHits_vs_trkChi2" : Hist(samp,cut,ele_trkHits,ele_chi2,storage=hist.storage.Weight()),
        #"sel_e2_numHits_vs_trkProb" : Hist(samp,cut,ele_trkHits,ele_prob,storage=hist.storage.Weight()),
        # 1D selected vertex histos
        "sel_vtx_type" : Hist(samp,cut,vtx_type,storage=hist.storage.Weight()),
        "sel_vtx_sign" : Hist(samp,cut,vtx_sign,storage=hist.storage.Weight()),
        "sel_vtx_dR" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_dR_over_pT" : Hist(samp,cut,dR_over_pT,storage=hist.storage.Weight()),
        "sel_vtx_dR_over_m" : Hist(samp,cut,dR_over_m,storage=hist.storage.Weight()),
        "sel_vtx_dR_over_pTm" : Hist(samp,cut,dR_over_pTm,storage=hist.storage.Weight()),
        "sel_vtx_dR_over_mpT" : Hist(samp,cut,dR_over_mpT,storage=hist.storage.Weight()),
        "sel_vtx_chi2" : Hist(samp,cut,vtx_chi2,storage=hist.storage.Weight()),
        "sel_vtx_prob" : Hist(samp,cut,vtx_prob,storage=hist.storage.Weight()),
        "sel_vtx_vxy" : Hist(samp,cut,vxy,storage=hist.storage.Weight()),
        "sel_vtx_vxy_zoom" : Hist(samp,cut,vxy_zoom,storage=hist.storage.Weight()),
        "sel_vtx_vxy_zoomzoom" : Hist(samp,cut,vxy_zoomzoom,storage=hist.storage.Weight()),
        "sel_vtx_vxySignif" : Hist(samp,cut,vtx_vxySignif,storage=hist.storage.Weight()),
        "sel_vtx_mass" : Hist(samp,cut,mass,storage=hist.storage.Weight()),
        "sel_vtx_minDxy" : Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        "sel_vtx_minDxy_fine": Hist(samp,cut,dxy_fine,storage=hist.storage.Weight()),
        "sel_vtx_pt" : Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_eta" : Hist(samp,cut,ele_eta,storage=hist.storage.Weight()),
        "sel_vtx_phi" : Hist(samp,cut,ele_phi,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso3" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso3M" : Hist(samp,cut,ele_PFRelIsoM,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso4" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso4M" : Hist(samp,cut,ele_PFRelIsoM,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso8" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFRelIso8M" : Hist(samp,cut,ele_PFRelIsoM,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso3" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso3M" : Hist(samp,cut,ele_PFIsoM,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso4" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso4M" : Hist(samp,cut,ele_PFIsoM,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso8" : Hist(samp,cut,ele_PFIso,storage=hist.storage.Weight()),
        #"sel_vtx_PFIso8M" : Hist(samp,cut,ele_PFIsoM,storage=hist.storage.Weight()),
        "sel_vtx_matchType" : Hist(samp,cut,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_max_chi2" : Hist(samp,cut,ele_chi2,storage=hist.storage.Weight()),
        "sel_vtx_min_pt" : Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_maxPFIso" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_minEledPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_matchType" : Hist(samp,cut,dphi,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj_vs_matchType" : Hist(samp,cut,dR,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_minEledPhiJ_vs_matchType" : Hist(samp,cut,dphi_generic,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "sel_vtx_mindRj_vs_matchType" : Hist(samp,cut,dR,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_mindPhiJ_vs_matchType" : Hist(samp,cut,dphi_generic,vtx_matchType,storage=hist.storage.Weight()),

        # 2D selected vertex histos
        "sel_vtx_mass_vs_mindxy" : Hist(samp,cut,mass,ele_dxy,storage=hist.storage.Weight()),
        "sel_vtx_mass_vs_vxy" : Hist(samp,cut,mass,vxy,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_chi2" : Hist(samp,cut,dphi,vtx_chi2,storage=hist.storage.Weight()),

        # Candidate ABCD planes
        "sel_vtx_METdPhi_vs_mindxy" : Hist(samp,cut,dphi,dxy,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_vxy" : Hist(samp,cut,dphi,vxy,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_mindPhiJ" : Hist(samp,cut,dphi,dphiJ,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_mindRj" : Hist(samp,cut,dphi,dRj,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_minEledPhiJ" : Hist(samp,cut,dphi,dphiJ,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_minEledRj" : Hist(samp,cut,dphi,dRj,storage=hist.storage.Weight()),
        "sel_vtx_maxPFIso_vs_mindxy" : Hist(samp,cut,pfiso,dxy,storage=hist.storage.Weight()),
        "sel_vtx_mindRj_vs_mindxy" : Hist(samp,cut,dR,dxy,storage=hist.storage.Weight()),
        "sel_vtx_mindRj_vs_maxPFIso" : Hist(samp,cut,dR,pfiso,storage=hist.storage.Weight()),
        "sel_vtx_mindPhiJ_vs_mindxy" : Hist(samp,cut,dphi,dxy,storage=hist.storage.Weight()),
        "sel_vtx_mindPhiJ_vs_maxPFIso" : Hist(samp,cut,dphi,pfiso,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj_vs_mindxy" : Hist(samp,cut,dR,dxy,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj_vs_maxPFIso" : Hist(samp,cut,dR,pfiso,storage=hist.storage.Weight()),
        "sel_vtx_minEledPhiJ_vs_mindxy" : Hist(samp,cut,dphi,dxy,storage=hist.storage.Weight()),
        "sel_vtx_minEledPhiJ_vs_maxPFIso" : Hist(samp,cut,dphi,pfiso,storage=hist.storage.Weight()),
        "sel_vtx_mindRj_vs_mindPhiJ" : Hist(samp,cut,dR,dphi,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj_vs_minEledPhiJ" : Hist(samp,cut,dR,dphi,storage=hist.storage.Weight()),
        
        # Misc other plots
        "met_over_lead_jet_pt" : Hist(samp,cut,met_over_pt,storage=hist.storage.Weight()),
        "lead_jet_met_dPhi" : Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "min_jet_met_dPhi" : Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "vtx_met_dPhi" : Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_over_pT" : Hist(samp,cut,METdPhi_over_pT,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_over_m" : Hist(samp,cut,METdPhi_over_m,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_over_pTm" : Hist(samp,cut,METdPhi_over_pTm,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_over_mpT" : Hist(samp,cut,METdPhi_over_mpT,storage=hist.storage.Weight()),
        "MET_pt" : Hist(samp,cut,met_pt,storage=hist.storage.Weight()),
        "nJets" : Hist(samp,cut,njets,storage=hist.storage.Weight()),
        "minBtag" : Hist(samp,cut,btag,storage=hist.storage.Weight()),
        "lead_jet_abseta" : Hist(samp,cut,jet_abseta,storage=hist.storage.Weight()),
        "lead_jet_pt" : Hist(samp,cut,jet_pt,storage=hist.storage.Weight()),
        #"dp_dotJet1" : Hist(samp,cut,angleDot,storage=hist.storage.Weight()),
        #"dp_dotJet12" : Hist(samp,cut,angleDot,storage=hist.storage.Weight()),
        # Gen plots
        "genEle_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genEle_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "genPos_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genPos_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "genEle_mindRjGen" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genEle_mindPhiJGen" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "genPos_mindRjGen" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genPos_mindPhiJGen" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "genEE_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genEE_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "genEE_mindRjGen" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "genEE_mindPhiJGen" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight())
    }
    return histograms

subroutines = []

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    e1 = events.sel_vtx.e1
    e2 = events.sel_vtx.e2
    min_dxy = np.minimum(np.abs(e1.dxy),np.abs(e2.dxy))
    max_pfiso = ak.where(e1.PFRelIso<e2.PFRelIso,e2.PFRelIso,e1.PFRelIso)
    wgt = events.eventWgt/sum_wgt
    vtx = events.sel_vtx
    # 1D selected electron 1 histos
    """histos["sel_e1_pt"].fill(samp=samp,cut=cut,pt=e1.pt,weight=wgt)
    histos["sel_e1_trkIso"].fill(samp=samp,cut=cut,trkIso=e1.trkIso,weight=wgt)
    histos["sel_e1_trkRelIso"].fill(samp=samp,cut=cut,relIso=e1.trkRelIso,weight=wgt)
    histos["sel_e1_PFRelIso3"].fill(samp=samp,cut=cut,relIso=e1.PFRelIso3,weight=wgt)
    histos["sel_e1_PFRelIso4"].fill(samp=samp,cut=cut,relIso=e1.PFRelIso4,weight=wgt)
    histos["sel_e1_PFRelIso8"].fill(samp=samp,cut=cut,relIso=e1.PFRelIso8,weight=wgt)
    histos["sel_e1_PFRelIso"].fill(samp=samp,cut=cut,relIso=e1.PFRelIso,weight=wgt)
    histos["sel_e1_PFIso3"].fill(samp=samp,cut=cut,iso=e1.PFIso3,weight=wgt)
    histos["sel_e1_PFIso4"].fill(samp=samp,cut=cut,iso=e1.PFIso4,weight=wgt)
    histos["sel_e1_PFIso8"].fill(samp=samp,cut=cut,iso=e1.PFIso8,weight=wgt)
    histos["sel_e1_trkChi2"].fill(samp=samp,cut=cut,chi2=e1.trkChi2,weight=wgt)
    histos["sel_e1_trkProb"].fill(samp=samp,cut=cut,prob=e1.trkProb,weight=wgt)
    histos["sel_e1_dxy"].fill(samp=samp,cut=cut,dxy=np.abs(e1.dxy),weight=wgt)
    histos["sel_e1_dxySignif"].fill(samp=samp,cut=cut,dxy_signif=np.abs(e1.dxy)/e1.dxyErr,weight=wgt)
    histos["sel_e1_angRes"].fill(samp=samp,cut=cut,angRes=e1.angRes,weight=wgt)
    histos["sel_e1_mindRj"].fill(samp=samp,cut=cut,dr=e1.mindRj,weight=wgt)
    histos["sel_e1_mindPhiJ"].fill(samp=samp,cut=cut,dphi=e1.mindPhiJ,weight=wgt)
    # 1D selected electron 2 histos
    histos["sel_e2_pt"].fill(samp=samp,cut=cut,pt=e2.pt,weight=wgt)
    histos["sel_e2_trkIso"].fill(samp=samp,cut=cut,trkIso=e2.trkIso,weight=wgt)
    histos["sel_e2_trkRelIso"].fill(samp=samp,cut=cut,relIso=e2.trkRelIso,weight=wgt)
    histos["sel_e2_PFRelIso3"].fill(samp=samp,cut=cut,relIso=e2.PFRelIso3,weight=wgt)
    histos["sel_e2_PFRelIso4"].fill(samp=samp,cut=cut,relIso=e2.PFRelIso4,weight=wgt)
    histos["sel_e2_PFRelIso8"].fill(samp=samp,cut=cut,relIso=e2.PFRelIso8,weight=wgt)
    histos["sel_e2_PFRelIso"].fill(samp=samp,cut=cut,relIso=e2.PFRelIso,weight=wgt)
    histos["sel_e2_PFIso3"].fill(samp=samp,cut=cut,iso=e2.PFIso3,weight=wgt)
    histos["sel_e2_PFIso4"].fill(samp=samp,cut=cut,iso=e2.PFIso4,weight=wgt)
    histos["sel_e2_PFIso8"].fill(samp=samp,cut=cut,iso=e2.PFIso8,weight=wgt)
    histos["sel_e2_trkChi2"].fill(samp=samp,cut=cut,chi2=e2.trkChi2,weight=wgt)
    histos["sel_e2_trkProb"].fill(samp=samp,cut=cut,prob=e2.trkProb,weight=wgt)
    histos["sel_e2_dxy"].fill(samp=samp,cut=cut,dxy=np.abs(e2.dxy),weight=wgt)
    histos["sel_e2_dxySignif"].fill(samp=samp,cut=cut,dxy_signif=np.abs(e2.dxy)/e2.dxyErr,weight=wgt)
    histos["sel_e2_angRes"].fill(samp=samp,cut=cut,angRes=e2.angRes,weight=wgt)
    histos["sel_e2_mindRj"].fill(samp=samp,cut=cut,dr=e2.mindRj,weight=wgt)
    histos["sel_e2_mindPhiJ"].fill(samp=samp,cut=cut,dphi=e2.mindPhiJ,weight=wgt)
    # 2D selected electron 1 histos
    histos["sel_e1_pt_vs_trkIso"].fill(samp=samp,cut=cut,pt=e1.pt,trkIso=e1.trkIso,weight=wgt)
    histos["sel_e1_pt_vs_trkRelIso"].fill(samp=samp,cut=cut,pt=e1.pt,relIso=e1.trkRelIso,weight=wgt)
    histos["sel_e1_pt_vs_chi2"].fill(samp=samp,cut=cut,pt=e1.pt,chi2=e1.trkChi2,weight=wgt)
    histos["sel_e1_chi2_vs_dxy"].fill(samp=samp,cut=cut,chi2=e1.trkChi2,dxy=np.abs(e1.dxy),weight=wgt)
    # 2D selected electron 2 histos
    histos["sel_e2_pt_vs_trkIso"].fill(samp=samp,cut=cut,pt=e2.pt,trkIso=e2.trkIso,weight=wgt)
    histos["sel_e2_pt_vs_trkRelIso"].fill(samp=samp,cut=cut,pt=e2.pt,relIso=e2.trkRelIso,weight=wgt)
    histos["sel_e2_pt_vs_chi2"].fill(samp=samp,cut=cut,pt=e2.pt,chi2=e2.trkChi2,weight=wgt)
    histos["sel_e2_chi2_vs_dxy"].fill(samp=samp,cut=cut,chi2=e2.trkChi2,dxy=np.abs(e2.dxy),weight=wgt)"""
    # 1D selected vertex histos
    histos["sel_vtx_type"].fill(samp=samp,cut=cut,type=vtx.typ,weight=wgt)
    histos["sel_vtx_sign"].fill(samp=samp,cut=cut,sign=vtx.sign,weight=wgt)
    histos["sel_vtx_dR"].fill(samp=samp,cut=cut,dr=vtx.dR,weight=wgt)
    histos["sel_vtx_dR_over_pT"].fill(samp=samp,cut=cut,dR_over_pT=vtx.dR/vtx.pt,weight=wgt)
    histos["sel_vtx_dR_over_m"].fill(samp=samp,cut=cut,dR_over_m=vtx.dR/vtx.m,weight=wgt)
    histos["sel_vtx_dR_over_pTm"].fill(samp=samp,cut=cut,dR_over_pTm=vtx.dR/(vtx.pt/vtx.m),weight=wgt)
    histos["sel_vtx_dR_over_mpT"].fill(samp=samp,cut=cut,dR_over_mpT=vtx.dR/(vtx.m/vtx.pt),weight=wgt)
    histos["sel_vtx_chi2"].fill(samp=samp,cut=cut,chi2=vtx.reduced_chi2,weight=wgt)
    histos["sel_vtx_prob"].fill(samp=samp,cut=cut,prob=vtx.prob,weight=wgt)
    histos["sel_vtx_vxy"].fill(samp=samp,cut=cut,vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_vxy_zoom"].fill(samp=samp,cut=cut,vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_vxy_zoomzoom"].fill(samp=samp,cut=cut,vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_vxySignif"].fill(samp=samp,cut=cut,vxy_signif=vtx.vxy/vtx.sigmavxy,weight=wgt)
    histos["sel_vtx_mass"].fill(samp=samp,cut=cut,mass=vtx.m,weight=wgt)
    histos["sel_vtx_minDxy"].fill(samp=samp,cut=cut,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_minDxy_fine"].fill(samp=samp,cut=cut,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_pt"].fill(samp=samp,cut=cut,pt=vtx.pt,weight=wgt)
    histos["sel_vtx_eta"].fill(samp=samp,cut=cut,eta=vtx.eta,weight=wgt)
    histos["sel_vtx_phi"].fill(samp=samp,cut=cut,phi=vtx.phi,weight=wgt)
    #histos["sel_vtx_PFRelIso3"].fill(samp=samp,cut=cut,relIso=vtx.PFRelIso3,weight=wgt)
    #histos["sel_vtx_PFRelIso3M"].fill(samp=samp,cut=cut,isoM=vtx.PFRelIso3*vtx.m,weight=wgt)
    #histos["sel_vtx_PFRelIso4"].fill(samp=samp,cut=cut,relIso=vtx.PFRelIso4,weight=wgt)
    #histos["sel_vtx_PFRelIso4M"].fill(samp=samp,cut=cut,isoM=vtx.PFRelIso4*vtx.m,weight=wgt)
    #histos["sel_vtx_PFRelIso8"].fill(samp=samp,cut=cut,relIso=vtx.PFRelIso8,weight=wgt)
    #histos["sel_vtx_PFRelIso8M"].fill(samp=samp,cut=cut,isoM=vtx.PFRelIso8*vtx.m,weight=wgt)
    #histos["sel_vtx_PFIso3"].fill(samp=samp,cut=cut,iso=vtx.PFIso3,weight=wgt)
    #histos["sel_vtx_PFIso3M"].fill(samp=samp,cut=cut,isoM=vtx.PFIso3*vtx.m,weight=wgt)
    #histos["sel_vtx_PFIso4"].fill(samp=samp,cut=cut,iso=vtx.PFIso4,weight=wgt)
    #histos["sel_vtx_PFIso4M"].fill(samp=samp,cut=cut,isoM=vtx.PFIso4*vtx.m,weight=wgt)
    #histos["sel_vtx_PFIso8"].fill(samp=samp,cut=cut,iso=vtx.PFIso8,weight=wgt)
    #histos["sel_vtx_PFIso8M"].fill(samp=samp,cut=cut,isoM=vtx.PFIso8*vtx.m,weight=wgt)
    histos["sel_vtx_max_chi2"].fill(samp=samp,cut=cut,chi2=ak.where(e1.trkChi2>e2.trkChi2,e1.trkChi2,e2.trkChi2),weight=wgt)
    histos["sel_vtx_min_pt"].fill(samp=samp,cut=cut,pt=ak.where(e1.pt<e2.pt,e1.pt,e2.pt),weight=wgt)
    histos["sel_vtx_maxPFIso"].fill(samp=samp,cut=cut,relIso=max_pfiso,weight=wgt)
    histos["sel_vtx_minEledRj"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),weight=wgt)
    histos["sel_vtx_minEledPhiJ"].fill(samp=samp,cut=cut,dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),weight=wgt)
    histos["sel_vtx_mindRj"].fill(samp=samp,cut=cut,dr=events.sel_vtx.mindRj,weight=wgt)
    histos["sel_vtx_mindPhiJ"].fill(samp=samp,cut=cut,dphi=events.sel_vtx.mindPhiJ,weight=wgt)
    # 2D selected vertex plots
    histos["sel_vtx_mass_vs_mindxy"].fill(samp=samp,cut=cut,mass=vtx.m,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_mass_vs_vxy"].fill(samp=samp,cut=cut,mass=vtx.m,vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_METdPhi_vs_chi2"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),chi2=vtx.reduced_chi2,weight=wgt)

    # Candidate ABCD planes
    histos["sel_vtx_METdPhi_vs_mindxy"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),dxy=min_dxy,weight=wgt)
    histos["sel_vtx_METdPhi_vs_vxy"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_METdPhi_vs_mindPhiJ"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),dphiJ=vtx.mindPhiJ,weight=wgt)
    histos["sel_vtx_METdPhi_vs_mindRj"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),drj=vtx.mindRj,weight=wgt)
    histos["sel_vtx_METdPhi_vs_minEledPhiJ"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),dphiJ=np.minimum(e1.mindPhiJ,e2.mindPhiJ),weight=wgt)
    histos["sel_vtx_METdPhi_vs_minEledRj"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),drj=np.minimum(e1.mindRj,e2.mindRj),weight=wgt)
    histos["sel_vtx_maxPFIso_vs_mindxy"].fill(samp=samp,cut=cut,pfiso=max_pfiso,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_mindRj_vs_mindxy"].fill(samp=samp,cut=cut,dr=vtx.mindRj,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_mindRj_vs_maxPFIso"].fill(samp=samp,cut=cut,dr=vtx.mindRj,pfiso=max_pfiso,weight=wgt)
    histos["sel_vtx_mindPhiJ_vs_mindxy"].fill(samp=samp,cut=cut,dphi=vtx.mindPhiJ,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_mindPhiJ_vs_maxPFIso"].fill(samp=samp,cut=cut,dphi=vtx.mindPhiJ,pfiso=max_pfiso,weight=wgt)
    histos["sel_vtx_minEledRj_vs_mindxy"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),dxy=min_dxy,weight=wgt)
    histos["sel_vtx_minEledRj_vs_maxPFIso"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),pfiso=max_pfiso,weight=wgt)
    histos["sel_vtx_minEledPhiJ_vs_mindxy"].fill(samp=samp,cut=cut,dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),dxy=min_dxy,weight=wgt)
    histos["sel_vtx_minEledPhiJ_vs_maxPFIso"].fill(samp=samp,cut=cut,dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),pfiso=max_pfiso,weight=wgt)
    histos["sel_vtx_mindRj_vs_mindPhiJ"].fill(samp=samp,cut=cut,dr=vtx.mindRj,dphi=vtx.mindPhiJ,weight=wgt)
    histos["sel_vtx_minEledRj_vs_minEledPhiJ"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),weight=wgt)

    # Misc other plots
    histos["met_over_lead_jet_pt"].fill(samp=samp,cut=cut,met_over_pt=events.PFMET.pt/events.PFJet.pt[:,0],weight=wgt)
    histos["lead_jet_met_dPhi"].fill(samp=samp,cut=cut,dphi=np.abs(events.PFJet.METdPhi[:,0]),weight=wgt)
    histos["min_jet_met_dPhi"].fill(samp=samp,cut=cut,dphi=ak.min(np.abs(events.PFJet.METdPhi),axis=1))
    histos["vtx_met_dPhi"].fill(samp=samp,cut=cut,dphi=np.abs(events.sel_vtx.METdPhi),weight=wgt)
    histos["sel_vtx_METdPhi_over_pT"].fill(samp=samp,cut=cut,METdPhi_over_pT=np.abs(events.sel_vtx.METdPhi)/vtx.pt,weight=wgt)
    histos["sel_vtx_METdPhi_over_m"].fill(samp=samp,cut=cut,METdPhi_over_m=np.abs(events.sel_vtx.METdPhi)/vtx.m,weight=wgt)
    histos["sel_vtx_METdPhi_over_pTm"].fill(samp=samp,cut=cut,METdPhi_over_pTm=np.abs(events.sel_vtx.METdPhi)/(vtx.pt/vtx.m),weight=wgt)
    histos["sel_vtx_METdPhi_over_mpT"].fill(samp=samp,cut=cut,METdPhi_over_mpT=np.abs(events.sel_vtx.METdPhi)/(vtx.m/vtx.pt),weight=wgt)
    histos["MET_pt"].fill(samp=samp,cut=cut,met_pt=events.PFMET.pt,weight=wgt)
    histos["nJets"].fill(samp=samp,cut=cut,njets=ak.count(events.PFJet.pt,axis=1),weight=wgt)
    histos["minBtag"].fill(samp=samp,cut=cut,btag=ak.fill_none(ak.min(events.PFJet.bTag,axis=1),-1),weight=wgt)
    histos["lead_jet_abseta"].fill(samp=samp,cut=cut,eta=np.abs(ak.fill_none(ak.pad_none(events.PFJet.eta,1),-999))[:,0],weight=wgt)
    histos["lead_jet_pt"].fill(samp=samp,cut=cut,pt=ak.fill_none(ak.pad_none(events.PFJet.pt,1),-999)[:,0],weight=wgt)
    #histos["dp_dotJet1"].fill(samp=samp,cut=cut,dot=events.DP_dotJet1,weight=wgt)
    #histos["dp_dotJet12"].fill(samp=samp,cut=cut,dot=events.DP_dotJet12,weight=wgt)

    if info["type"] == "signal":
        histos["sel_vtx_matchType"].fill(samp=samp,cut=cut,mtype=vtx.match,weight=wgt)
        histos["sel_vtx_minEledRj_vs_matchType"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),mtype=vtx.match,weight=wgt)
        histos["sel_vtx_minEledPhiJ_vs_matchType"].fill(samp=samp,cut=cut,dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),mtype=vtx.match,weight=wgt)
        histos["sel_vtx_mindRj_vs_matchType"].fill(samp=samp,cut=cut,dr=events.sel_vtx.mindRj,mtype=vtx.match,weight=wgt)
        histos["sel_vtx_mindPhiJ_vs_matchType"].fill(samp=samp,cut=cut,dphi=events.sel_vtx.mindPhiJ,mtype=vtx.match,weight=wgt)
        histos["sel_vtx_METdPhi_vs_matchType"].fill(samp=samp,cut=cut,dphi=np.abs(events.sel_vtx.METdPhi),mtype=vtx.match,weight=wgt)
        
        histos["genEle_mindRj"].fill(samp=samp,cut=cut,dr=events.GenEle.mindRj,weight=wgt)
        histos["genEle_mindPhiJ"].fill(samp=samp,cut=cut,dphi=events.GenEle.mindPhiJ,weight=wgt)
        histos["genEle_mindRjGen"].fill(samp=samp,cut=cut,dr=events.GenEle.mindRjGen,weight=wgt)
        histos["genEle_mindPhiJGen"].fill(samp=samp,cut=cut,dphi=events.GenEle.mindPhiJGen,weight=wgt)

        histos["genPos_mindRj"].fill(samp=samp,cut=cut,dr=events.GenPos.mindRj,weight=wgt)
        histos["genPos_mindPhiJ"].fill(samp=samp,cut=cut,dphi=events.GenPos.mindPhiJ,weight=wgt)
        histos["genPos_mindRjGen"].fill(samp=samp,cut=cut,dr=events.GenPos.mindRjGen,weight=wgt)
        histos["genPos_mindPhiJGen"].fill(samp=samp,cut=cut,dphi=events.GenPos.mindPhiJGen,weight=wgt)

        histos["genEE_mindRj"].fill(samp=samp,cut=cut,dr=events.genEE.mindRj,weight=wgt)
        histos["genEE_mindPhiJ"].fill(samp=samp,cut=cut,dphi=events.genEE.mindPhiJ,weight=wgt)
        histos["genEE_mindRjGen"].fill(samp=samp,cut=cut,dr=events.genEE.mindRjGen,weight=wgt)
        histos["genEE_mindPhiJGen"].fill(samp=samp,cut=cut,dphi=events.genEE.mindPhiJGen,weight=wgt)

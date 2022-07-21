from histobins import *
from coffea.hist import Hist
import numpy as np
import awkward as ak

histograms = {
    # 1D selected electron histos
    "sel_e_pt" : Hist("Events",sample,cut,ele_pt),
    "sel_e_trkIso" : Hist("Events",sample,cut,ele_trkIso),
    "sel_e_trkRelIso" : Hist("Events",sample,cut,ele_trkRelIso),
    "sel_e_trkChi2" : Hist("Events",sample,cut,ele_chi2),
    "sel_e_trkProb" : Hist("Events",sample,cut,ele_prob),
    "sel_e_dxy" : Hist("Events",sample,cut,ele_dxy),
    "sel_e_dxySignif" : Hist("Events",sample,cut,ele_dxySignif),
    "sel_e_angRes" : Hist("Events",sample,cut,ele_angRes),
    # 2D selected electron histos
    "sel_e_pt_vs_trkIso" : Hist("Events",sample,cut,ele_pt,ele_trkIso),
    "sel_e_pt_vs_trkRelIso" : Hist("Events",sample,cut,ele_pt,ele_trkRelIso),
    "sel_e_pt_vs_chi2" : Hist("Events",sample,cut,ele_pt,ele_chi2),
    "sel_e_chi2_vs_dxy" : Hist("Events",sample,cut,ele_chi2,ele_dxy),
    # 1D selected vertex histos
    "sel_vtx_dR" : Hist("Events",sample,cut,dR),
    "sel_vtx_chi2" : Hist("Events",sample,cut,vtx_chi2),
    "sel_vtx_prob" : Hist("Events",sample,cut,vtx_prob),
    "sel_vtx_vxy" : Hist("Events",sample,cut,vxy),
    "sel_vtx_vxy_zoom" : Hist("Events",sample,cut,vxy_zoom),
    "sel_vtx_vxy_zoomzoom" : Hist("Events",sample,cut,vxy_zoomzoom),
    "sel_vtx_vxySignif" : Hist("Events",sample,cut,vtx_vxySignif),
    "sel_vtx_mass" : Hist("Events",sample,cut,mass),
    "sel_vtx_minDxy" : Hist("Events",sample,cut,ele_dxy),
    "sel_vtx_pt" : Hist("Events",sample,cut,ele_pt),
    # 2D selected vertex histos
    "sel_vtx_pt_vs_vxy" : Hist("Events",sample,cut,vtx_pt,vxy),
    "sel_vtx_pt_vs_vxy_zoom" : Hist("Events",sample,cut,vtx_pt,vxy_zoom),
    "sel_vtx_pt_vs_vxy_zoomzoom" : Hist("Events",sample,cut,vtx_pt,vxy_zoomzoom),
    "sel_vtx_dR_vs_vxy" : Hist("Events",sample,cut,dR,vxy),
    "sel_vtx_dR_vs_vxy_zoom" : Hist("Events",sample,cut,dR,vxy_zoom),
    "sel_vtx_dR_vs_vxy_zoomzoom" : Hist("Events",sample,cut,dR,vxy_zoomzoom),
    "sel_vtx_chi2_vs_vxy" : Hist("Events",sample,cut,vtx_chi2,vxy),
    "sel_vtx_chi2_vs_vxy_zoom" : Hist("Events",sample,cut,vtx_chi2,vxy_zoom),
    "sel_vtx_chi2_vs_vxy_zoomzoom" : Hist("Events",sample,cut,vtx_chi2,vxy_zoomzoom),
    # Misc other plots
    "jet_met_dPhi" : Hist("Events",sample,cut,jet_met_dPhi)
}

subroutines = []

def fillHistos(events,histos,samp,cut):
    sel_e = ak.concatenate((events.sel_e1,events.sel_e2))
    min_dxy = ak.min(ak.concatenate((events.sel_e1.dxy[:,np.newaxis],events.sel_e2.dxy[:,np.newaxis]),axis=1),axis=1)
    # 1D selected electron histos
    histos["sel_e_pt"].fill(sample=samp,cut=cut,pt=sel_e.pt)
    histos["sel_e_trkIso"].fill(sample=samp,cut=cut,trkIso=sel_e.trkIso)
    histos["sel_e_trkRelIso"].fill(sample=samp,cut=cut,relIso=sel_e.trkRelIso)
    histos["sel_e_trkChi2"].fill(sample=samp,cut=cut,chi2=sel_e.trkChi2)
    histos["sel_e_trkProb"].fill(sample=samp,cut=cut,prob=sel_e.trkProb)
    histos["sel_e_dxy"].fill(sample=samp,cut=cut,dxy=sel_e.dxy)
    histos["sel_e_dxySignif"].fill(sample=samp,cut=cut,dxy_signif=sel_e.dxy/sel_e.dxyErr)
    histos["sel_e_angRes"].fill(sample=samp,cut=cut,angRes=sel_e.angRes)
    # 2D selected electron histos
    histos["sel_e_pt_vs_trkIso"].fill(sample=samp,cut=cut,pt=sel_e.pt,trkIso=sel_e.trkIso)
    histos["sel_e_pt_vs_trkRelIso"].fill(sample=samp,cut=cut,pt=sel_e.pt,relIso=sel_e.trkRelIso)
    histos["sel_e_pt_vs_chi2"].fill(sample=samp,cut=cut,pt=sel_e.pt,chi2=sel_e.trkChi2)
    histos["sel_e_chi2_vs_dxy"].fill(sample=samp,cut=cut,chi2=sel_e.trkChi2,dxy=sel_e.dxy)
    # 1D selected vertex histos
    histos["sel_vtx_dR"].fill(sample=samp,cut=cut,dR=events.sel_vtx.dR)
    histos["sel_vtx_chi2"].fill(sample=samp,cut=cut,chi2=events.sel_vtx.reduced_chi2)
    histos["sel_vtx_prob"].fill(sample=samp,cut=cut,prob=events.sel_vtx.prob)
    histos["sel_vtx_vxy"].fill(sample=samp,cut=cut,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_vxy_zoom"].fill(sample=samp,cut=cut,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_vxy_zoomzoom"].fill(sample=samp,cut=cut,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_vxySignif"].fill(sample=samp,cut=cut,vxy_signif=events.sel_vtx.vxy/events.sel_vtx.sigmavxy)
    histos["sel_vtx_mass"].fill(sample=samp,cut=cut,mass=events.sel_vtx.m)
    histos["sel_vtx_minDxy"].fill(sample=samp,cut=cut,dxy=min_dxy)
    histos["sel_vtx_pt"].fill(sample=samp,cut=cut,pt=events.sel_vtx.pt)
    # 2D selected vertex histos
    histos["sel_vtx_pt_vs_vxy"].fill(sample=samp,cut=cut,pt=events.sel_vtx.pt,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_pt_vs_vxy_zoom"].fill(sample=samp,cut=cut,pt=events.sel_vtx.pt,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_pt_vs_vxy_zoomzoom"].fill(sample=samp,cut=cut,pt=events.sel_vtx.pt,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_dR_vs_vxy"].fill(sample=samp,cut=cut,dR=events.sel_vtx.dR,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_dR_vs_vxy_zoom"].fill(sample=samp,cut=cut,dR=events.sel_vtx.dR,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_dR_vs_vxy_zoomzoom"].fill(sample=samp,cut=cut,dR=events.sel_vtx.dR,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_chi2_vs_vxy"].fill(sample=samp,cut=cut,chi2=events.sel_vtx.reduced_chi2,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_chi2_vs_vxy_zoom"].fill(sample=samp,cut=cut,chi2=events.sel_vtx.reduced_chi2,vxy=events.sel_vtx.vxy)
    histos["sel_vtx_chi2_vs_vxy_zoomzoom"].fill(sample=samp,cut=cut,chi2=events.sel_vtx.reduced_chi2,vxy=events.sel_vtx.vxy)
    # Misc other plots
    histos["jet_met_dPhi"].fill(sample=samp,cut=cut,jet_met_dphi=ak.flatten(events.JetMETdPhi))

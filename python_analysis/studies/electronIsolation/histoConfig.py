import sys
sys.path.append("../../configs/histo_configs/")
from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        "sel_vtx_type" : Hist(samp,cut,vtx_type,storage=hist.storage.Weight()),
        "sel_vtx_sign" : Hist(samp,cut,vtx_sign,storage=hist.storage.Weight()),
        "sel_vtx_dR" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
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
        "sel_vtx_matchType" : Hist(samp,cut,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_max_chi2" : Hist(samp,cut,ele_chi2,storage=hist.storage.Weight()),
        "sel_vtx_min_pt" : Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_maxPFIso" : Hist(samp,cut,ele_PFRelIso,storage=hist.storage.Weight()),
        "sel_vtx_minEledRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_minEledPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight()),
        "sel_vtx_mindRj" : Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_mindPhiJ" : Hist(samp,cut,dphi_generic,storage=hist.storage.Weight())
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
    # 1D selected vertex histos
    histos["sel_vtx_type"].fill(samp=samp,cut=cut,type=vtx.typ,weight=wgt)
    histos["sel_vtx_sign"].fill(samp=samp,cut=cut,sign=vtx.sign,weight=wgt)
    histos["sel_vtx_dR"].fill(samp=samp,cut=cut,dr=vtx.dR,weight=wgt)
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
    histos["sel_vtx_max_chi2"].fill(samp=samp,cut=cut,chi2=ak.where(e1.trkChi2>e2.trkChi2,e1.trkChi2,e2.trkChi2),weight=wgt)
    histos["sel_vtx_min_pt"].fill(samp=samp,cut=cut,pt=ak.where(e1.pt<e2.pt,e1.pt,e2.pt),weight=wgt)
    histos["sel_vtx_maxPFIso"].fill(samp=samp,cut=cut,relIso=max_pfiso,weight=wgt)
    histos["sel_vtx_minEledRj"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),weight=wgt)
    histos["sel_vtx_minEledPhiJ"].fill(samp=samp,cut=cut,dphi=np.minimum(e1.mindPhiJ,e2.mindPhiJ),weight=wgt)
    histos["sel_vtx_mindRj"].fill(samp=samp,cut=cut,dr=events.sel_vtx.mindRj,weight=wgt)
    histos["sel_vtx_mindPhiJ"].fill(samp=samp,cut=cut,dphi=events.sel_vtx.mindPhiJ,weight=wgt)

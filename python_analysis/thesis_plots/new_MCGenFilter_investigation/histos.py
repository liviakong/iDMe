from histobins import *
from hist import Hist
from hist.axis import StrCategory, Regular, Integer, IntCategory
import hist
import numpy as np
import awkward as ak

# General Purpose
samp = StrCategory([],name="samp",label="Sample Name",growth=True)
cut = StrCategory([],name="cut",label="Cut Applied",growth=True)

# functions to make histograms
def parse_axis(a):
    name = a[0]
    if type(a[1]) == list:
        assert len(a) == 2
        if type(a[1][0]) == str:
            axis = StrCategory(a[1],name=name,label=name)
        else:
            axis = IntCategory(a[1],name=name,label=name)
    else:
        assert len(a) == 4
        axis = Regular(a[1],a[2],a[3],name=name,label=name)
    return axis

def histo(name,*args):
    axes = [samp,cut]
    for ax in args:
        axes.append(parse_axis(ax))
    return Hist(*axes,storage=hist.storage.Weight())
    

def make_histograms():
    kind = ('kind',['L','R'])
    match = ('match',[0,1])
    histograms = {
        # quantities associated w/ gen objects
        "gen_met" : histo("gen_met",('met',100,50,300)),
        "gen_met_noWgt" : histo("gen_met",('met',100,50,300)),
        "gen_dR" : histo("gen_dR",('dR',250,0,5)),
        "gen_vxy" : histo("gen_vxy",('vxy',500,0,50)),
        "gen_leadpT" : histo("gen_leadpT",("pt",100,0,50)),
        "gen_vtxMETdPhi" : histo("gen_vtxMETdPhi",("dphi",64,-3.2,3.2)),
        "gen_jetMETdPhi" : histo("gen_jetMETdPhi",("dphi",64,-3.2,3.2)),
        "gen_ele_pt" : histo('gen_ele_pt',('pt',100,0,50)),
        #
        "match_ele_pt" : histo("match_ele_pt",kind,('pt',100,0,50)),
        "match_ele_dxy" : histo("match_ele_dxy",kind,('dxy',50,0,5)),
        "match_ele_dxySignif" : histo("match_ele_dxySignif",kind,('signif',100,0,20)),
        "match_ele_trkChi2" : histo("match_ele_trkChi2",kind,('chi2',100,0,5)),
        "match_ele_trkRelIso" : histo("match_ele_trkRelIso",kind,('iso',100,0,5)),
        "match_ele_calRelIso" : histo("match_ele_calRelIso",kind,('iso',100,0,5)),
        "match_ele_PFRelIso" : histo("match_ele_PFRelIso",kind,('iso',100,0,5)),
        "match_ele_mindRJets" : histo("match_ele_mindRJets",kind,('dR',100,0,5)),
        "match_ele_sigmaIetaIeta" : histo("match_ele_sigmaIetaIeta",kind,('sieie',100,0,0.1)),
        "match_ele_absdEtaSeed" : histo("match_ele_absdEtaSeed",kind,('detaSeed',100,0,5)),
        "match_ele_absdPhiIn" : histo("match_ele_absdPhiIn",kind,('dphi',100,0,3)),
        "match_ele_HOverE" : histo("match_ele_HOverE",kind,('hoe',100,0,200)),
        "match_ele_1E1p" : histo("match_ele_1E1p",kind,('emp',100,0,5)),
        "match_ele_expMissing" : histo("match_ele_expMissing",kind,('missing',10,0,10)),
        "match_ele_passConvVeto" : histo("match_ele_passConvVeto",kind,('passVeto',[0,1])),
        #
        "fake_ele_pt" : histo("fake_ele_pt",kind,('pt',100,0,50)),
        "fake_ele_dxy" : histo("fake_ele_dxy",kind,('dxy',50,0,5)),
        "fake_ele_dxySignif" : histo("fake_ele_dxySignif",kind,('signif',100,0,20)),
        "fake_ele_trkChi2" : histo("fake_ele_trkChi2",kind,('chi2',100,0,5)),
        "fake_ele_trkRelIso" : histo("fake_ele_trkRelIso",kind,('iso',100,0,5)),
        "fake_ele_calRelIso" : histo("fake_ele_calRelIso",kind,('iso',100,0,5)),
        "fake_ele_PFRelIso" : histo("fake_ele_PFRelIso",kind,('iso',100,0,5)),
        "fake_ele_mindRJets" : histo("fake_ele_mindRJets",kind,('dR',100,0,5)),
        "fake_ele_sigmaIetaIeta" : histo("fake_ele_sigmaIetaIeta",kind,('sieie',100,0,0.1)),
        "fake_ele_absdEtaSeed" : histo("fake_ele_absdEtaSeed",kind,('detaSeed',100,0,5)),
        "fake_ele_absdPhiIn" : histo("fake_ele_absdPhiIn",kind,('dphi',100,0,3)),
        "fake_ele_HOverE" : histo("fake_ele_HOverE",kind,('hoe',100,0,200)),
        "fake_ele_1E1p" : histo("fake_ele_1E1p",kind,('emp',100,0,5)),
        "fake_ele_expMissing" : histo("fake_ele_expMissing",kind,('missing',10,0,10)),
        "fake_ele_passConvVeto" : histo("fake_ele_passConvVeto",kind,('passVeto',[0,1])),
        #
        "match_vtx_dR" : histo("match_vtx_dR",('dR',250,0,5)),
        "match_vtx_mindxy" : histo("match_vtx_mindxy",('dxy',50,0,5)),
        "match_vtx_vxy" : histo("match_vtx_vxy",('vxy',500,0,50)),
        "match_vtx_vxy_relDiff" : histo("match_vtx_vxy_relDiff",('vxy_diff',50,-2,2)),
        "match_vtx_leadpT" : histo("match_vtx_leadpT",("pt",50,0,50)),
        "match_vtx_METdPhi" : histo("match_vtx_METdPhi",("dphi",64,-3.2,3.2)),
        #"match_vtx_mindRj" : histo("match_vtx_mindRj",('dR',250,0,5)),
        #
        "fake_vtx_dR" : histo("fake_vtx_dR",('dR',250,0,5)),
        "fake_vtx_mindxy" : histo("fake_vtx_mindxy",('dxy',50,0,5)),
        "fake_vtx_vxy" : histo("fake_vtx_vxy",('vxy',500,0,50)),
        "fake_vtx_leadpT" : histo("fake_vtx_leadpT",("pt",50,0,50)),
        "fake_vtx_METdPhi" : histo("fake_vtx_METdPhi",("dphi",64,-3.2,3.2)),
        #"fake_vtx_mindRj" : histo("fake_vtx_mindRj",('dR',250,0,5)),
         #                    
        "PFMET" : histo("PFMET",("met",100,200,350))
    }
    return histograms

subroutines = []

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    wgt = events.eventWgt/sum_wgt
    
    if info["type"] == "signal":
        match_pf = ak.flatten(events.Electron[events.Electron.genMatched])
        fake_pf = ak.flatten(events.Electron[~events.Electron.genMatched])
        match_lpt = ak.flatten(events.LptElectron[events.LptElectron.genMatched])
        fake_lpt = ak.flatten(events.LptElectron[~events.LptElectron.genMatched])
        match_vtx = ak.flatten(events.vtx[events.vtx.isMatched])
        fake_vtx = ak.flatten(events.vtx[events.vtx.isMatched])
        hasMatchVtx = ak.any(events.vtx.isMatched,axis=1)

        histos["gen_met_noWgt"].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=1)
        histos['gen_met'].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=wgt)
        histos["gen_dR"].fill(samp=samp,cut=cut,dR=events.genEE.dr,weight=wgt)
        histos["gen_vxy"].fill(samp=samp,cut=cut,vxy=events.GenEle.vxy,weight=wgt)
        histos["gen_leadpT"].fill(samp=samp,cut=cut,pt=np.maximum(events.GenEle.pt,events.GenPos.pt),weight=wgt)
        histos["gen_vtxMETdPhi"].fill(samp=samp,cut=cut,dphi=events.genEE.METdPhi,weight=wgt)
        histos["gen_jetMETdPhi"].fill(samp=samp,cut=cut,dphi=events.GenJetMETdPhi,weight=wgt)
        histos["gen_ele_pt"].fill(samp=samp,cut=cut,pt=events.GenEle.pt,weight=wgt)
        histos["gen_ele_pt"].fill(samp=samp,cut=cut,pt=events.GenPos.pt,weight=wgt)
        #
        histos["match_ele_pt"].fill(samp=samp,cut=cut,kind='R',pt=match_pf.pt,weight=1)
        histos["match_ele_dxy"].fill(samp=samp,cut=cut,kind='R',dxy=match_pf.dxy,weight=1)
        histos["match_ele_dxySignif"].fill(samp=samp,cut=cut,kind='R',signif=match_pf.dxy/match_pf.dxyErr,weight=1)
        histos["match_ele_trkChi2"].fill(samp=samp,cut=cut,kind='R',chi2=match_pf.trkChi2,weight=1)
        histos["match_ele_trkRelIso"].fill(samp=samp,cut=cut,kind='R',iso=match_pf.trkRelIso,weight=1)
        histos["match_ele_calRelIso"].fill(samp=samp,cut=cut,kind='R',iso=match_pf.calRelIso,weight=1)
        histos["match_ele_PFRelIso"].fill(samp=samp,cut=cut,kind='R',iso=match_pf.PFRelIso,weight=1)
        histos["match_ele_mindRJets"].fill(samp=samp,cut=cut,kind='R',dR=ak.min(match_pf.dRJets,axis=1),weight=1)
        histos["match_ele_sigmaIetaIeta"].fill(samp=samp,cut=cut,kind='R',sieie=match_pf.full55sigmaIetaIeta,weight=1)
        histos["match_ele_absdEtaSeed"].fill(samp=samp,cut=cut,kind='R',detaSeed=match_pf.absdEtaSeed,weight=1)
        histos["match_ele_absdPhiIn"].fill(samp=samp,cut=cut,kind='R',dphi=match_pf.absdPhiIn,weight=1)
        histos["match_ele_HOverE"].fill(samp=samp,cut=cut,kind='R',hoe=match_pf.HoverE,weight=1)
        histos["match_ele_1E1p"].fill(samp=samp,cut=cut,kind='R',emp=match_pf.abs1overEm1overP,weight=1)
        histos["match_ele_expMissing"].fill(samp=samp,cut=cut,kind='R',missing=match_pf.expMissingInnerHits,weight=1)
        histos["match_ele_passConvVeto"].fill(samp=samp,cut=cut,kind='R',passVeto=match_pf.conversionVeto,weight=1)

        histos["match_ele_pt"].fill(samp=samp,cut=cut,kind='L',pt=match_lpt.pt,weight=1)
        histos["match_ele_dxy"].fill(samp=samp,cut=cut,kind='L',dxy=match_lpt.dxy,weight=1)
        histos["match_ele_dxySignif"].fill(samp=samp,cut=cut,kind='L',signif=match_lpt.dxy/match_lpt.dxyErr,weight=1)
        histos["match_ele_trkChi2"].fill(samp=samp,cut=cut,kind='L',chi2=match_lpt.trkChi2,weight=1)
        histos["match_ele_trkRelIso"].fill(samp=samp,cut=cut,kind='L',iso=match_lpt.trkRelIso,weight=1)
        histos["match_ele_calRelIso"].fill(samp=samp,cut=cut,kind='L',iso=match_lpt.calRelIso,weight=1)
        histos["match_ele_PFRelIso"].fill(samp=samp,cut=cut,kind='L',iso=match_lpt.PFRelIso,weight=1)
        histos["match_ele_mindRJets"].fill(samp=samp,cut=cut,kind='L',dR=ak.min(match_lpt.dRJets,axis=1),weight=1)
        histos["match_ele_sigmaIetaIeta"].fill(samp=samp,cut=cut,kind='L',sieie=match_lpt.full55sigmaIetaIeta,weight=1)
        histos["match_ele_absdEtaSeed"].fill(samp=samp,cut=cut,kind='L',detaSeed=match_lpt.absdEtaSeed,weight=1)
        histos["match_ele_absdPhiIn"].fill(samp=samp,cut=cut,kind='L',dphi=match_lpt.absdPhiIn,weight=1)
        histos["match_ele_HOverE"].fill(samp=samp,cut=cut,kind='L',hoe=match_lpt.HoverE,weight=1)
        histos["match_ele_1E1p"].fill(samp=samp,cut=cut,kind='L',emp=match_lpt.abs1overEm1overP,weight=1)
        histos["match_ele_expMissing"].fill(samp=samp,cut=cut,kind='L',missing=match_lpt.expMissingInnerHits,weight=1)
        histos["match_ele_passConvVeto"].fill(samp=samp,cut=cut,kind='L',passVeto=match_lpt.conversionVeto,weight=1)
        #
        histos["fake_ele_pt"].fill(samp=samp,cut=cut,kind='R',pt=fake_pf.pt,weight=1)
        histos["fake_ele_dxy"].fill(samp=samp,cut=cut,kind='R',dxy=fake_pf.dxy,weight=1)
        histos["fake_ele_dxySignif"].fill(samp=samp,cut=cut,kind='R',signif=fake_pf.dxy/fake_pf.dxyErr,weight=1)
        histos["fake_ele_trkChi2"].fill(samp=samp,cut=cut,kind='R',chi2=fake_pf.trkChi2,weight=1)
        histos["fake_ele_trkRelIso"].fill(samp=samp,cut=cut,kind='R',iso=fake_pf.trkRelIso,weight=1)
        histos["fake_ele_calRelIso"].fill(samp=samp,cut=cut,kind='R',iso=fake_pf.calRelIso,weight=1)
        histos["fake_ele_PFRelIso"].fill(samp=samp,cut=cut,kind='R',iso=fake_pf.PFRelIso,weight=1)
        histos["fake_ele_mindRJets"].fill(samp=samp,cut=cut,kind='R',dR=ak.min(fake_pf.dRJets,axis=1),weight=1)
        histos["fake_ele_sigmaIetaIeta"].fill(samp=samp,cut=cut,kind='R',sieie=fake_pf.full55sigmaIetaIeta,weight=1)
        histos["fake_ele_absdEtaSeed"].fill(samp=samp,cut=cut,kind='R',detaSeed=fake_pf.absdEtaSeed,weight=1)
        histos["fake_ele_absdPhiIn"].fill(samp=samp,cut=cut,kind='R',dphi=fake_pf.absdPhiIn,weight=1)
        histos["fake_ele_HOverE"].fill(samp=samp,cut=cut,kind='R',hoe=fake_pf.HoverE,weight=1)
        histos["fake_ele_1E1p"].fill(samp=samp,cut=cut,kind='R',emp=fake_pf.abs1overEm1overP,weight=1)
        histos["fake_ele_expMissing"].fill(samp=samp,cut=cut,kind='R',missing=fake_pf.expMissingInnerHits,weight=1)
        histos["fake_ele_passConvVeto"].fill(samp=samp,cut=cut,kind='R',passVeto=fake_pf.conversionVeto,weight=1)

        histos["fake_ele_pt"].fill(samp=samp,cut=cut,kind='L',pt=fake_lpt.pt,weight=1)
        histos["fake_ele_dxy"].fill(samp=samp,cut=cut,kind='L',dxy=fake_lpt.dxy,weight=1)
        histos["fake_ele_dxySignif"].fill(samp=samp,cut=cut,kind='L',signif=fake_lpt.dxy/fake_lpt.dxyErr,weight=1)
        histos["fake_ele_trkChi2"].fill(samp=samp,cut=cut,kind='L',chi2=fake_lpt.trkChi2,weight=1)
        histos["fake_ele_trkRelIso"].fill(samp=samp,cut=cut,kind='L',iso=fake_lpt.trkRelIso,weight=1)
        histos["fake_ele_calRelIso"].fill(samp=samp,cut=cut,kind='L',iso=fake_lpt.calRelIso,weight=1)
        histos["fake_ele_PFRelIso"].fill(samp=samp,cut=cut,kind='L',iso=fake_lpt.PFRelIso,weight=1)
        histos["fake_ele_mindRJets"].fill(samp=samp,cut=cut,kind='L',dR=ak.min(fake_lpt.dRJets,axis=1),weight=1)
        histos["fake_ele_sigmaIetaIeta"].fill(samp=samp,cut=cut,kind='L',sieie=fake_lpt.full55sigmaIetaIeta,weight=1)
        histos["fake_ele_absdEtaSeed"].fill(samp=samp,cut=cut,kind='L',detaSeed=fake_lpt.absdEtaSeed,weight=1)
        histos["fake_ele_absdPhiIn"].fill(samp=samp,cut=cut,kind='L',dphi=fake_lpt.absdPhiIn,weight=1)
        histos["fake_ele_HOverE"].fill(samp=samp,cut=cut,kind='L',hoe=fake_lpt.HoverE,weight=1)
        histos["fake_ele_1E1p"].fill(samp=samp,cut=cut,kind='L',emp=fake_lpt.abs1overEm1overP,weight=1)
        histos["fake_ele_expMissing"].fill(samp=samp,cut=cut,kind='L',missing=fake_lpt.expMissingInnerHits,weight=1)
        histos["fake_ele_passConvVeto"].fill(samp=samp,cut=cut,kind='L',passVeto=fake_lpt.conversionVeto,weight=1)
        #
        histos["match_vtx_dR"].fill(samp=samp,cut=cut,dR=match_vtx.dR,weight=1)
        histos["match_vtx_mindxy"].fill(samp=samp,cut=cut,dxy=np.minimum(match_vtx.e1.dxy,match_vtx.e2.dxy),weight=1)
        histos["match_vtx_vxy"].fill(samp=samp,cut=cut,vxy=match_vtx.vxy,weight=1)
        histos["match_vtx_vxy_relDiff"].fill(samp=samp,cut=cut,vxy_diff=(match_vtx.vxy-events.GenEle.vxy[hasMatchVtx])/events.GenEle.vxy[hasMatchVtx],weight=1)
        histos["match_vtx_leadpT"].fill(samp=samp,cut=cut,pt=np.maximum(match_vtx.e1.pt,match_vtx.e2.pt),weight=1)
        histos["match_vtx_METdPhi"].fill(samp=samp,cut=cut,dphi=match_vtx.METdPhi,weight=1)
        #histos["match_vtx_mindRj"].fill(samp=samp,cut=cut,dR=ak.min(match_vtx.dRJets,axis=1),weight=1)

        histos["fake_vtx_dR"].fill(samp=samp,cut=cut,dR=fake_vtx.dR,weight=1)
        histos["fake_vtx_mindxy"].fill(samp=samp,cut=cut,dxy=np.minimum(fake_vtx.e1.dxy,fake_vtx.e2.dxy),weight=1)
        histos["fake_vtx_vxy"].fill(samp=samp,cut=cut,vxy=fake_vtx.vxy,weight=1)
        histos["fake_vtx_leadpT"].fill(samp=samp,cut=cut,pt=np.maximum(fake_vtx.e1.pt,fake_vtx.e2.pt),weight=1)
        histos["fake_vtx_METdPhi"].fill(samp=samp,cut=cut,dphi=fake_vtx.METdPhi,weight=1)
        #histos["fake_vtx_mindRj"].fill(samp=samp,cut=cut,dR=ak.min(fake_vtx.dRJets,axis=1),weight=1)
         #                    
        histos["PFMET"].fill(samp=samp,cut=cut,met=events.PFMET.pt,weight=wgt)
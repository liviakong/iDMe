from hist import Hist
from hist.axis import StrCategory, Regular, Integer, IntCategory
import hist
import numpy as np
import awkward as ak

# General Purpose
samp = StrCategory([],name="samp",label="Sample Name",growth=True)
cut = StrCategory([],name="cut",label="Cut Applied",growth=True)

# functions to make histograms
class myHisto:
    def __init__(self):
        self.histograms = {}
        self.samp = "NO_SAMPLE"
        self.cut = "NO_CUT"

        # axis compendium
        self.samp = StrCategory([],name="samp",label="Sample Name",growth=True)
        self.cut = StrCategory([],name="cut",label="Cut Applied",growth=True)
        self.ele_type = self.parse_axis(('ele_type',['L','R']))
        self.match_type = self.parse_axis(('match_type',['L','R']))
        self.match = self.parse_axis(('match',[0,1]))
        self.met = self.parse_axis(('met',100,50,300))
        self.dR = self.parse_axis(('dR',100,0,5))
        self.vxy1 = self.parse_axis(('vxy',100,0,1))
        self.vxy10 = self.parse_axis(('vxy',100,0,10))
        self.vxy100 = self.parse_axis(('vxy',100,0,100))
        self.ele_pt = self.parse_axis(("pt",50,0,50))
        self.dphi = self.parse_axis(("phi",64,-3.2,3.2))
        self.phi = self.parse_axis(("phi",64,-3.2,3.2))
        self.abs_dphi = self.parse_axis(('phi',100,0,3))
        self.eta = self.parse_axis(('eta',100,-2.5,2.5))
        self.dxy = self.parse_axis(('dxy',500,0,5))
        self.r3 = self.parse_axis(('r3',500,0,100))
        self.dxy_signif = self.parse_axis(('signif',100,0,20))
        self.trk_chi2 = self.parse_axis(('chi2',100,0,5))
        self.vtx_chi2 = self.parse_axis(('chi2',100,0,5))
        self.iso = self.parse_axis(('iso',100,0,5))
        self.sieie = self.parse_axis(('sieie',100,0,0.1))
        self.detaseed = self.parse_axis(('detaSeed',100,0,5))
        self.HoverE = self.parse_axis(('hoe',100,0,200))
        self.EinvMinusPinv = self.parse_axis(('emp',100,0,5))
        self.numMissingHits = self.parse_axis(('missing',10,0,10))
        self.passConvVeto = self.parse_axis(('passVeto',[0,1]))
        self.vxy_relDiff = self.parse_axis(('rel_diff',50,-2,2))
        self.isMatched = self.parse_axis(('matched',[0,1]))
        self.isPF = self.parse_axis(('isPF',[0,1]))
        self.numHits = self.parse_axis(('numHits',20,0,20))
        self.trkProb = self.parse_axis(('prob',100,0,1))
        self.IDScore = self.parse_axis(('id',100,-1,3))
        self.ele_passID = self.parse_axis(('passID',[0,1,2,3]))
        self.ele_IDtype = self.parse_axis(("IDtype",['Base','dR']))
        self.vtx_type = self.parse_axis(('vtype',['LL','LR','RR']))
        self.vtx_mass = self.parse_axis(('mass',100,0,30))
        self.vtx_sign = self.parse_axis(('sign',[-1,1]))
        self.vtx_pt = self.parse_axis(('pt',100,0,50))
        self.ele_ptRes = self.parse_axis(('ptres',100,-2,2))
        self.sigReco = self.parse_axis(('reco',[0,1]))
        self.vtxMatch = self.parse_axis(('match',[0,1]))
        self.dRCategories = self.parse_axis(('dRCat',['0to0p1','0p1to0p5','0p5toInf']))
        self.vxyCategories = self.parse_axis(('vxyCat',['0to1','1to5','5to10','10to15','15toInf']))
        self.ptCategories = self.parse_axis(('ptCat',['0to5','5to10','10to20','20toInf']))

    def make(self,name,*args,**hist_kwargs):
        if name in self.histograms.keys():
            print(f"Histogram {name} already exists! Skipping")
            return
        axes = [self.samp,self.cut]
        for ax in args:
            if type(ax) == tuple:
                axes.append(self.parse_axis(ax))
            else:
                axes.append(getattr(self,ax))
        self.histograms[name] = Hist(*axes,storage=hist.storage.Weight(),**hist_kwargs)
    
    def fill(self,name,**kwargs):
        self.histograms[name].fill(samp=self.samp,cut=self.cut,**kwargs)
    
    def parse_axis(self,a):
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

def make_histograms():
    h = myHisto()
    
    # quantities associated w/ gen objects
    h.make("gen_met",'met')
    h.make("gen_dR",'dR')
    h.make("gen_vxy1",'vxy1')
    h.make("gen_vxy10",'vxy10')
    h.make("gen_vxy100",'vxy100')
    h.make("gen_leadpT",'ele_pt')
    h.make("gen_vtx_METdPhi",'dphi')
    h.make("gen_jetMETdPhi",'dphi')
    h.make("gen_vtx_mass",'vtx_mass')
    h.make('gen_vtx_pt','vtx_pt')
    h.make('gen_vtx_eta','eta')
    h.make('gen_vtx_phi','phi')

    h.make("gen_met_unwgt",'met')
    h.make("gen_dR_unwgt",'dR')
    h.make("gen_vxy1_unwgt",'vxy1')
    h.make("gen_vxy10_unwgt",'vxy10')
    h.make("gen_vxy100_unwgt",'vxy100')
    h.make("gen_leadpT_unwgt",'ele_pt')
    h.make("gen_vtx_METdPhi_unwgt",'dphi')
    h.make("gen_jetMETdPhi_unwgt",'dphi')
    h.make("gen_vtx_mass_unwgt",'vtx_mass')
    h.make('gen_vtx_pt_unwgt','vtx_pt')
    h.make('gen_vtx_eta_unwgt','eta')
    h.make('gen_vtx_phi_unwgt','phi')

    h.make('gen_ele_pt','ele_pt')
    h.make('gen_ele_eta','eta')
    h.make('gen_ele_phi','phi')
    h.make('gen_ele_vxy1','vxy1')
    h.make('gen_ele_vxy10','vxy10')
    h.make('gen_ele_vxy100','vxy100')
    h.make('gen_ele_pt_dRbin_vxybin','ele_pt','dRCategories','vxyCategories')
    h.make('gen_ele_dR_ptbin_vxybin','dR','ptCategories','vxyCategories')
    h.make('gen_ele_vxy1_dRbin_ptbin','vxy1','dRCategories','ptCategories')
    h.make('gen_ele_vxy10_dRbin_ptbin','vxy10','dRCategories','ptCategories')
    h.make('gen_ele_vxy100_dRbin_ptbin','vxy100','dRCategories','ptCategories')

    # reco --> vertexing efficiency
    h.make("signalReco_vs_vtxMatch",'sigReco','vtxMatch')
    h.make("signalReco_vs_vtxMatch_unwgt",'sigReco','vtxMatch')
    
    # matched reco electron, reco variables
    h.make("match_ele_pt",'match_type','ele_passID','ele_IDtype','ele_pt')
    h.make("match_ele_eta",'match_type','ele_passID','ele_IDtype','eta')
    h.make("match_ele_phi",'match_type','ele_passID','ele_IDtype','phi')
    h.make("match_ele_dxy",'match_type','ele_passID','ele_IDtype','dxy')
    h.make("match_ele_dxySignif",'match_type','ele_passID','ele_IDtype','dxy_signif')
    h.make("match_ele_trkChi2",'match_type','ele_passID','ele_IDtype','trk_chi2')
    h.make('match_ele_trkProb','match_type','ele_passID','ele_IDtype','trkProb')
    h.make("match_ele_trkRelIso",'match_type','ele_passID','ele_IDtype','iso')
    h.make("match_ele_calRelIso",'match_type','ele_passID','ele_IDtype','iso')
    h.make("match_ele_PFRelIso",'match_type','ele_passID','ele_IDtype','iso')
    h.make("match_ele_miniRelIso",'match_type','ele_passID','ele_IDtype','iso')
    h.make("match_ele_mindRJets",'match_type','ele_passID','ele_IDtype','dR')
    h.make("match_ele_isPF",'match_type','ele_passID','ele_IDtype','isPF')
    h.make("match_ele_numTrkHits",'match_type','ele_passID','ele_IDtype','numHits')
    h.make("match_ele_numPixHits",'match_type','ele_passID','ele_IDtype','numHits')
    h.make("match_ele_numStripHits",'match_type','ele_passID','ele_IDtype','numHits')
    h.make("match_ele_mindPhiJets",'match_type','ele_passID','ele_IDtype','dphi')
    h.make("match_ele_sigmaIetaIeta",'match_type','ele_passID','ele_IDtype','sieie')
    h.make("match_ele_absdEtaSeed",'match_type','ele_passID','ele_IDtype','detaseed')
    h.make("match_ele_absdPhiIn",'match_type','ele_passID','ele_IDtype','abs_dphi')
    h.make("match_ele_HOverE",'match_type','ele_passID','ele_IDtype','HoverE')
    h.make("match_ele_1E1p",'match_type','ele_passID','ele_IDtype','EinvMinusPinv')
    h.make("match_ele_expMissing",'match_type','ele_passID','ele_IDtype','numMissingHits')
    h.make("match_ele_passConvVeto",'match_type','ele_passID','ele_IDtype','passConvVeto')
    h.make("match_ele_IDscore",'match_type','ele_passID','ele_IDtype','IDScore')
    
    # matched reco electrons, corresponding gen object variables
    h.make("match_ele_gen_pt",'match_type','ele_passID','ele_IDtype','ele_pt')
    h.make('match_ele_gen_vxy1','match_type','ele_passID','ele_IDtype','vxy1')
    h.make('match_ele_gen_vxy10','match_type','ele_passID','ele_IDtype','vxy10')
    h.make('match_ele_gen_vxy100','match_type','ele_passID','ele_IDtype','vxy100')
    h.make('match_ele_gen_pt_dRbin_vxybin','match_type','ele_passID','ele_IDtype','ele_pt','dRCategories','vxyCategories')
    h.make('match_ele_gen_dR_ptbin_vxybin','match_type','ele_passID','ele_IDtype','dR','ptCategories','vxyCategories')
    h.make('match_ele_gen_vxy1_dRbin_ptbin','match_type','ele_passID','ele_IDtype','vxy1','dRCategories','ptCategories')
    h.make('match_ele_gen_vxy10_dRbin_ptbin','match_type','ele_passID','ele_IDtype','vxy10','dRCategories','ptCategories')
    h.make('match_ele_gen_vxy100_dRbin_ptbin','match_type','ele_passID','ele_IDtype','vxy100','dRCategories','ptCategories')
    h.make("match_ele_gen_eta",'match_type','ele_passID','ele_IDtype','eta')
    h.make('match_ele_gen_phi','match_type','ele_passID','ele_IDtype','phi')
    h.make('match_ele_gen_dR','match_type','ele_passID','ele_IDtype','dR')
    h.make('match_ele_ptRes','match_type','ele_passID','ele_IDtype','ele_ptRes')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    wgt = events.eventWgt/sum_wgt
    
    if info["type"] == "signal":
        # defining stuff
        hasMatch_pf = (ak.count(events.Electron.pt,axis=1)>0) &\
                      (ak.count_nonzero(events.Electron.genMatched,axis=1)>0)

        hasMatch_lpt = (ak.count(events.LptElectron.pt,axis=1)>0) &\
                       (ak.count_nonzero(events.LptElectron.genMatched,axis=1)>0)

        match_pf = events[hasMatch_pf].Electron
        match_pf = match_pf[match_pf.genMatched]
        genObj_pf = ak.where(match_pf.matchType==-1,events[hasMatch_pf].GenEle,events[hasMatch_pf].GenPos)
        match_pf = ak.flatten(match_pf)
        genObj_pf = ak.flatten(genObj_pf)
        match_pf_passID = ak.values_astype(match_pf.passID,int)
        match_pf_passIDBasic = ak.values_astype(match_pf.passIDBasic,int)

        match_lpt = events[hasMatch_lpt].LptElectron
        match_lpt = match_lpt[match_lpt.genMatched]
        genObj_lpt = ak.where(match_lpt.matchType==-1,events[hasMatch_lpt].GenEle,events[hasMatch_lpt].GenPos)
        match_lpt = ak.flatten(match_lpt)
        genObj_lpt = ak.flatten(genObj_lpt)
        match_lpt_passID = ak.values_astype(match_lpt.passID,int)
        match_lpt_passIDBasic = ak.values_astype(match_lpt.passIDBasic,int)

        fake_pf = ak.flatten(events.Electron[~events.Electron.genMatched])
        fake_lpt = ak.flatten(events.LptElectron[~events.LptElectron.genMatched])

        fake_pf_passID = ak.values_astype(fake_pf.passID,int)
        fake_lpt_passID = ak.values_astype(fake_lpt.passID,int)

        match_vtx = ak.flatten(events.vtx[events.vtx.isMatched])
        match_vtx_genObj = events.genEE[ak.any(events.vtx.isMatched,axis=1)]
        fake_vtx = ak.flatten(events.vtx[~events.vtx.isMatched])

        hasMatchVtx = ak.values_astype(ak.any(events.vtx.isMatched,axis=1),int)
        signalReconstructed = ak.values_astype(events.signalReconstructed,int)

        ### FILLING HISTOGRAMS ###

        h.fill("gen_met",met=events.GenMET.pt,weight=wgt)
        h.fill("gen_dR",dR=events.genEE.dr,weight=wgt)
        h.fill("gen_vxy1",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_vxy10",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_vxy100",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_leadpT",pt=np.maximum(events.GenEle.pt,events.GenPos.pt),weight=wgt)
        h.fill("gen_vtx_METdPhi",phi=events.genEE.METdPhi,weight=wgt)
        h.fill("gen_jetMETdPhi",phi=events.GenJetMETdPhi,weight=wgt)
        h.fill("gen_vtx_mass",mass=events.genEE.mass,weight=wgt)
        h.fill('gen_vtx_pt',pt=events.genEE.pt,weight=wgt)
        h.fill('gen_vtx_eta',eta=events.genEE.eta,weight=wgt)
        h.fill('gen_vtx_phi',phi=events.genEE.phi,weight=wgt)
        
        
        h.fill("gen_met_unwgt",met=events.GenMET.pt,weight=1)
        h.fill("gen_dR_unwgt",dR=events.genEE.dr,weight=1)
        h.fill("gen_vxy1_unwgt",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_vxy10_unwgt",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_vxy100_unwgt",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_leadpT_unwgt",pt=np.maximum(events.GenEle.pt,events.GenPos.pt),weight=1)
        h.fill("gen_vtx_METdPhi_unwgt",phi=events.genEE.METdPhi,weight=1)
        h.fill("gen_jetMETdPhi_unwgt",phi=events.GenJetMETdPhi,weight=1)
        h.fill("gen_vtx_mass_unwgt",mass=events.genEE.mass,weight=1)
        h.fill('gen_vtx_pt_unwgt',pt=events.genEE.pt,weight=1)
        h.fill('gen_vtx_eta_unwgt',eta=events.genEE.eta,weight=1)
        h.fill('gen_vtx_phi_unwgt',phi=events.genEE.phi,weight=1)
        
        h.fill("gen_ele_pt",pt=events.GenEle.pt,weight=1)
        h.fill("gen_ele_pt",pt=events.GenPos.pt,weight=1)
        h.fill("gen_ele_vxy1",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy1",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_vxy10",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy10",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_vxy100",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy100",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_eta",eta=events.GenEle.eta,weight=1)
        h.fill("gen_ele_eta",eta=events.GenPos.eta,weight=1)
        h.fill("gen_ele_phi",phi=events.GenEle.phi,weight=1)
        h.fill("gen_ele_phi",phi=events.GenPos.phi,weight=1)
        h.fill('gen_ele_pt_dRbin_vxybin',pt=events.GenEle.pt,dRCat=events.GenEle.dRbin,vxyCat=events.GenEle.vxyBin,weight=1)
        h.fill('gen_ele_dR_ptbin_vxybin',dR=events.GenEle.dr,ptCat=events.GenEle.ptBin,vxyCat=events.GenEle.vxyBin,weight=1)
        h.fill('gen_ele_vxy1_dRbin_ptbin',vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,ptCat=events.GenEle.ptBin,weight=1)
        h.fill('gen_ele_vxy10_dRbin_ptbin',vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,ptCat=events.GenEle.ptBin,weight=1)
        h.fill('gen_ele_vxy100_dRbin_ptbin',vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,ptCat=events.GenEle.ptBin,weight=1)
        h.fill('gen_ele_pt_dRbin_vxybin',pt=events.GenPos.pt,dRCat=events.GenPos.dRbin,vxyCat=events.GenPos.vxyBin,weight=1)
        h.fill('gen_ele_dR_ptbin_vxybin',dR=events.GenPos.dr,ptCat=events.GenPos.ptBin,vxyCat=events.GenPos.vxyBin,weight=1)
        h.fill('gen_ele_vxy1_dRbin_ptbin',vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,ptCat=events.GenPos.ptBin,weight=1)
        h.fill('gen_ele_vxy10_dRbin_ptbin',vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,ptCat=events.GenPos.ptBin,weight=1)
        h.fill('gen_ele_vxy100_dRbin_ptbin',vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,ptCat=events.GenPos.ptBin,weight=1)
        
        #
        h.fill('signalReco_vs_vtxMatch',reco=signalReconstructed,match=hasMatchVtx,weight=wgt)
        h.fill('signalReco_vs_vtxMatch_unwgt',reco=signalReconstructed,match=hasMatchVtx,weight=1)  
        
        #
        h.fill("match_ele_pt",IDtype='dR',match_type='R',passID=match_pf_passID,pt=match_pf.pt,weight=1)
        h.fill("match_ele_eta",IDtype='dR',match_type='R',passID=match_pf_passID,eta=match_pf.eta,weight=1)
        h.fill("match_ele_phi",IDtype='dR',match_type='R',passID=match_pf_passID,phi=match_pf.phi,weight=1)
        h.fill("match_ele_dxy",IDtype='dR',match_type='R',passID=match_pf_passID,dxy=match_pf.dxy,weight=1)
        h.fill("match_ele_dxySignif",IDtype='dR',match_type='R',passID=match_pf_passID,signif=match_pf.dxy/match_pf.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",IDtype='dR',match_type='R',passID=match_pf_passID,chi2=match_pf.trkChi2,weight=1)
        h.fill("match_ele_trkProb",IDtype='dR',match_type='R',passID=match_pf_passID,prob=match_pf.trkProb)
        h.fill("match_ele_trkRelIso",IDtype='dR',match_type='R',passID=match_pf_passID,iso=match_pf.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",IDtype='dR',match_type='R',passID=match_pf_passID,iso=match_pf.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",IDtype='dR',match_type='R',passID=match_pf_passID,iso=match_pf.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",IDtype='dR',match_type='R',passID=match_pf_passID,iso=match_pf.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",IDtype='dR',match_type='R',passID=match_pf_passID,dR=match_pf.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",IDtype='dR',match_type='R',passID=match_pf_passID,phi=match_pf.mindPhiJ,weight=1)
        h.fill('match_ele_isPF',IDtype='dR',match_type='R',passID=match_pf_passID,isPF=ak.values_astype(match_pf.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',IDtype='dR',match_type='R',passID=match_pf_passID,numHits=match_pf.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',IDtype='dR',match_type='R',passID=match_pf_passID,numHits=match_pf.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',IDtype='dR',match_type='R',passID=match_pf_passID,numHits=match_pf.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",IDtype='dR',match_type='R',passID=match_pf_passID,sieie=match_pf.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",IDtype='dR',match_type='R',passID=match_pf_passID,detaSeed=match_pf.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",IDtype='dR',match_type='R',passID=match_pf_passID,phi=match_pf.absdPhiIn,weight=1)
        h.fill("match_ele_HOverE",IDtype='dR',match_type='R',passID=match_pf_passID,hoe=match_pf.HoverE,weight=1)
        h.fill("match_ele_1E1p",IDtype='dR',match_type='R',passID=match_pf_passID,emp=match_pf.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",IDtype='dR',match_type='R',passID=match_pf_passID,missing=match_pf.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",IDtype='dR',match_type='R',passID=match_pf_passID,passVeto=match_pf.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",IDtype='dR',match_type='R',passID=match_pf_passID,id=match_pf.IDscore,weight=1)

        h.fill("match_ele_pt",IDtype='dR',match_type='L',passID=match_lpt_passID,pt=match_lpt.pt,weight=1)
        h.fill("match_ele_eta",IDtype='dR',match_type='L',passID=match_lpt_passID,eta=match_lpt.eta,weight=1)
        h.fill("match_ele_phi",IDtype='dR',match_type='L',passID=match_lpt_passID,phi=match_lpt.phi,weight=1)
        h.fill("match_ele_dxy",IDtype='dR',match_type='L',passID=match_lpt_passID,dxy=match_lpt.dxy,weight=1)
        h.fill("match_ele_dxySignif",IDtype='dR',match_type='L',passID=match_lpt_passID,signif=match_lpt.dxy/match_lpt.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",IDtype='dR',match_type='L',passID=match_lpt_passID,chi2=match_lpt.trkChi2,weight=1)
        h.fill("match_ele_trkProb",IDtype='dR',match_type='L',passID=match_lpt_passID,prob=match_lpt.trkProb)
        h.fill("match_ele_trkRelIso",IDtype='dR',match_type='L',passID=match_lpt_passID,iso=match_lpt.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",IDtype='dR',match_type='L',passID=match_lpt_passID,iso=match_lpt.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",IDtype='dR',match_type='L',passID=match_lpt_passID,iso=match_lpt.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",IDtype='dR',match_type='L',passID=match_lpt_passID,iso=match_lpt.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",IDtype='dR',match_type='L',passID=match_lpt_passID,dR=match_lpt.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",IDtype='dR',match_type='L',passID=match_lpt_passID,phi=match_lpt.mindPhiJ,weight=1)
        h.fill('match_ele_isPF',IDtype='dR',match_type='L',passID=match_lpt_passID,isPF=ak.values_astype(match_lpt.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',IDtype='dR',match_type='L',passID=match_lpt_passID,numHits=match_lpt.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',IDtype='dR',match_type='L',passID=match_lpt_passID,numHits=match_lpt.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',IDtype='dR',match_type='L',passID=match_lpt_passID,numHits=match_lpt.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",IDtype='dR',match_type='L',passID=match_lpt_passID,sieie=match_lpt.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",IDtype='dR',match_type='L',passID=match_lpt_passID,detaSeed=match_lpt.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",IDtype='dR',match_type='L',passID=match_lpt_passID,phi=match_lpt.absdPhiIn,weight=1)
        h.fill("match_ele_HOverE",IDtype='dR',match_type='L',passID=match_lpt_passID,hoe=match_lpt.HoverE,weight=1)
        h.fill("match_ele_1E1p",IDtype='dR',match_type='L',passID=match_lpt_passID,emp=match_lpt.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",IDtype='dR',match_type='L',passID=match_lpt_passID,missing=match_lpt.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",IDtype='dR',match_type='L',passID=match_lpt_passID,passVeto=match_lpt.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",IDtype='dR',match_type='L',passID=match_lpt_passID,id=match_lpt.IDscore,weight=1)
        
        #
        h.fill("match_ele_gen_pt",IDtype='dR',match_type='R',passID=match_pf_passID,pt=genObj_pf.pt,weight=1)
        h.fill("match_ele_gen_vxy1",IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,weight=1)
        h.fill("match_ele_gen_vxy10",IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,weight=1)
        h.fill("match_ele_gen_vxy100",IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,weight=1)
        h.fill('match_ele_gen_pt_dRbin_vxybin',IDtype='dR',match_type='R',passID=match_pf_passID,pt=genObj_pf.pt,dRCat=genObj_pf.dRbin,vxyCat=genObj_pf.vxyBin)
        h.fill('match_ele_gen_dR_ptbin_vxybin',IDtype='dR',match_type='R',passID=match_pf_passID,dR=genObj_pf.dr,ptCat=genObj_pf.ptBin,vxyCat=genObj_pf.vxyBin)
        h.fill('match_ele_gen_vxy1_dRbin_ptbin',IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill('match_ele_gen_vxy10_dRbin_ptbin',IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill('match_ele_gen_vxy100_dRbin_ptbin',IDtype='dR',match_type='R',passID=match_pf_passID,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill("match_ele_gen_eta",IDtype='dR',match_type='R',passID=match_pf_passID,eta=genObj_pf.eta,weight=1)
        h.fill("match_ele_gen_phi",IDtype='dR',match_type='R',passID=match_pf_passID,phi=genObj_pf.phi,weight=1)
        h.fill("match_ele_gen_dR",IDtype='dR',match_type='R',passID=match_pf_passID,dR=genObj_pf.dr,weight=1)
        h.fill("match_ele_ptRes",IDtype='dR',match_type='R',passID=match_pf_passID,ptres=(genObj_pf.pt-match_pf.pt)/genObj_pf.pt,weight=1)

        h.fill("match_ele_gen_pt",IDtype='dR',match_type='L',passID=match_lpt_passID,pt=genObj_lpt.pt,weight=1)
        h.fill("match_ele_gen_vxy1",IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,weight=1)
        h.fill("match_ele_gen_vxy10",IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,weight=1)
        h.fill("match_ele_gen_vxy100",IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,weight=1)
        h.fill('match_ele_gen_pt_dRbin_vxybin',IDtype='dR',match_type='L',passID=match_lpt_passID,pt=genObj_lpt.pt,dRCat=genObj_lpt.dRbin,vxyCat=genObj_lpt.vxyBin)
        h.fill('match_ele_gen_dR_ptbin_vxybin',IDtype='dR',match_type='L',passID=match_lpt_passID,dR=genObj_lpt.dr,ptCat=genObj_lpt.ptBin,vxyCat=genObj_lpt.vxyBin)
        h.fill('match_ele_gen_vxy1_dRbin_ptbin',IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill('match_ele_gen_vxy10_dRbin_ptbin',IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill('match_ele_gen_vxy100_dRbin_ptbin',IDtype='dR',match_type='L',passID=match_lpt_passID,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill("match_ele_gen_eta",IDtype='dR',match_type='L',passID=match_lpt_passID,eta=genObj_lpt.eta,weight=1)
        h.fill("match_ele_gen_phi",IDtype='dR',match_type='L',passID=match_lpt_passID,phi=genObj_lpt.phi,weight=1)
        h.fill("match_ele_gen_dR",IDtype='dR',match_type='L',passID=match_lpt_passID,dR=genObj_lpt.dr,weight=1)
        h.fill("match_ele_ptRes",IDtype='dR',match_type='L',passID=match_lpt_passID,ptres=(genObj_lpt.pt-match_lpt.pt)/genObj_lpt.pt,weight=1)

        ################################################
        ### with basic id (no dR(e,jet) requirement) ###
        ################################################
        h.fill("match_ele_pt",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,pt=match_pf.pt,weight=1)
        h.fill("match_ele_eta",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,eta=match_pf.eta,weight=1)
        h.fill("match_ele_phi",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,phi=match_pf.phi,weight=1)
        h.fill("match_ele_dxy",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,dxy=match_pf.dxy,weight=1)
        h.fill("match_ele_dxySignif",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,signif=match_pf.dxy/match_pf.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,chi2=match_pf.trkChi2,weight=1)
        h.fill("match_ele_trkProb",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,prob=match_pf.trkProb)
        h.fill("match_ele_trkRelIso",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,iso=match_pf.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,iso=match_pf.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,iso=match_pf.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,iso=match_pf.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,dR=match_pf.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,phi=match_pf.mindPhiJ,weight=1)
        h.fill('match_ele_isPF',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,isPF=ak.values_astype(match_pf.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,numHits=match_pf.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,numHits=match_pf.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,numHits=match_pf.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,sieie=match_pf.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,detaSeed=match_pf.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,phi=match_pf.absdPhiIn,weight=1)
        h.fill("match_ele_HOverE",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,hoe=match_pf.HoverE,weight=1)
        h.fill("match_ele_1E1p",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,emp=match_pf.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,missing=match_pf.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,passVeto=match_pf.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,id=match_pf.IDscore,weight=1)

        h.fill("match_ele_pt",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,pt=match_lpt.pt,weight=1)
        h.fill("match_ele_eta",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,eta=match_lpt.eta,weight=1)
        h.fill("match_ele_phi",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,phi=match_lpt.phi,weight=1)
        h.fill("match_ele_dxy",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,dxy=match_lpt.dxy,weight=1)
        h.fill("match_ele_dxySignif",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,signif=match_lpt.dxy/match_lpt.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,chi2=match_lpt.trkChi2,weight=1)
        h.fill("match_ele_trkProb",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,prob=match_lpt.trkProb)
        h.fill("match_ele_trkRelIso",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,iso=match_lpt.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,iso=match_lpt.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,iso=match_lpt.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,iso=match_lpt.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,dR=match_lpt.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,phi=match_lpt.mindPhiJ,weight=1)
        h.fill('match_ele_isPF',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,isPF=ak.values_astype(match_lpt.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,numHits=match_lpt.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,numHits=match_lpt.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,numHits=match_lpt.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,sieie=match_lpt.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,detaSeed=match_lpt.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,phi=match_lpt.absdPhiIn,weight=1)
        h.fill("match_ele_HOverE",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,hoe=match_lpt.HoverE,weight=1)
        h.fill("match_ele_1E1p",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,emp=match_lpt.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,missing=match_lpt.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,passVeto=match_lpt.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,id=match_lpt.IDscore,weight=1)
        
        #
        h.fill("match_ele_gen_pt",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,pt=genObj_pf.pt,weight=1)
        h.fill("match_ele_gen_vxy1",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,weight=1)
        h.fill("match_ele_gen_vxy10",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,weight=1)
        h.fill("match_ele_gen_vxy100",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,weight=1)
        h.fill('match_ele_gen_pt_dRbin_vxybin',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,pt=genObj_pf.pt,dRCat=genObj_pf.dRbin,vxyCat=genObj_pf.vxyBin)
        h.fill('match_ele_gen_dR_ptbin_vxybin',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,dR=genObj_pf.dr,ptCat=genObj_pf.ptBin,vxyCat=genObj_pf.vxyBin)
        h.fill('match_ele_gen_vxy1_dRbin_ptbin',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill('match_ele_gen_vxy10_dRbin_ptbin',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill('match_ele_gen_vxy100_dRbin_ptbin',IDtype='Base',match_type='R',passID=match_pf_passIDBasic,vxy=genObj_pf.vxy,dRCat=genObj_pf.dRbin,ptCat=genObj_pf.ptBin)
        h.fill("match_ele_gen_eta",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,eta=genObj_pf.eta,weight=1)
        h.fill("match_ele_gen_phi",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,phi=genObj_pf.phi,weight=1)
        h.fill("match_ele_gen_dR",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,dR=genObj_pf.dr,weight=1)
        h.fill("match_ele_ptRes",IDtype='Base',match_type='R',passID=match_pf_passIDBasic,ptres=(genObj_pf.pt-match_pf.pt)/genObj_pf.pt,weight=1)

        h.fill("match_ele_gen_pt",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,pt=genObj_lpt.pt,weight=1)
        h.fill("match_ele_gen_vxy1",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,weight=1)
        h.fill("match_ele_gen_vxy10",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,weight=1)
        h.fill("match_ele_gen_vxy100",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,weight=1)
        h.fill('match_ele_gen_pt_dRbin_vxybin',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,pt=genObj_lpt.pt,dRCat=genObj_lpt.dRbin,vxyCat=genObj_lpt.vxyBin)
        h.fill('match_ele_gen_dR_ptbin_vxybin',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,dR=genObj_lpt.dr,ptCat=genObj_lpt.ptBin,vxyCat=genObj_lpt.vxyBin)
        h.fill('match_ele_gen_vxy1_dRbin_ptbin',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill('match_ele_gen_vxy10_dRbin_ptbin',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill('match_ele_gen_vxy100_dRbin_ptbin',IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,vxy=genObj_lpt.vxy,dRCat=genObj_lpt.dRbin,ptCat=genObj_lpt.ptBin)
        h.fill("match_ele_gen_eta",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,eta=genObj_lpt.eta,weight=1)
        h.fill("match_ele_gen_phi",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,phi=genObj_lpt.phi,weight=1)
        h.fill("match_ele_gen_dR",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,dR=genObj_lpt.dr,weight=1)
        h.fill("match_ele_ptRes",IDtype='Base',match_type='L',passID=match_lpt_passIDBasic,ptres=(genObj_lpt.pt-match_lpt.pt)/genObj_lpt.pt,weight=1)
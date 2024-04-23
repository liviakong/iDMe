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
        self.ele_type = self.parse_axis(('ele_type',['L','R','Both']))
        self.match_type = self.parse_axis(('match_type',['L','R','Both']))
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
        self.sigReco = self.parse_axis(('sigReco',[0,1]))
        self.sigRecoPassID = self.parse_axis(('sigRecoPassID',[0,1]))
        self.vtxMatch = self.parse_axis(('vtxMatch',[0,1]))
        self.vtxPassID = self.parse_axis(('vtxPassID',[0,1]))
        self.dRCategories = self.parse_axis(('dRCat',['0to0p1','0p1to0p5','0p5toInf']))
        self.vxyCategories = self.parse_axis(('vxyCat',['0to1','1to5','5to10','10to15','15toInf']))
        self.ptCategories = self.parse_axis(('ptCat',['0to5','5to10','10to20','20toInf']))
        self.genEmatched = self.parse_axis(('eMatch',[0,1]))
        self.genPmatched = self.parse_axis(('pMatch',[0,1]))
        self.genEpassID = self.parse_axis(('eID',[0,1]))
        self.genPpassID = self.parse_axis(('pID',[0,1]))
        

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
    
    # reco --> vertexing efficiency
    h.make("signalReco_vs_vtxMatch_BaseID",'sigReco','vtxMatch','sigRecoPassID','vtxPassID')
    h.make("signalReco_vs_vtxMatch_dRID",'sigReco','vtxMatch','sigRecoPassID','vtxPassID')
    h.make("genEPMatched_BaseID",'genEmatched','genPmatched','genEpassID','genPpassID')
    h.make("genEPMatched_dRID",'genEmatched','genPmatched','genEpassID','genPpassID')
    h.make("genEPMatched_vtxMatched_BaseID",'genEmatched','genPmatched','genEpassID','genPpassID','vtxMatch','vtxPassID')
    h.make("genEPMatched_vtxMatched_dRID",'genEmatched','genPmatched','genEpassID','genPpassID','vtxMatch','vtxPassID')
    
    h.make("signalReco_vs_vtxMatch_BaseID_unwgt",'sigReco','vtxMatch','sigRecoPassID','vtxPassID')
    h.make("signalReco_vs_vtxMatch_dRID_unwgt",'sigReco','vtxMatch','sigRecoPassID','vtxPassID')
    h.make("genEPMatched_BaseID_unwgt",'genEmatched','genPmatched','genEpassID','genPpassID')
    h.make("genEPMatched_dRID_unwgt",'genEmatched','genPmatched','genEpassID','genPpassID')
    h.make("genEPMatched_vtxMatched_BaseID_unwgt",'genEmatched','genPmatched','genEpassID','genPpassID','vtxMatch','vtxPassID')
    h.make("genEPMatched_vtxMatched_dRID_unwgt",'genEmatched','genPmatched','genEpassID','genPpassID','vtxMatch','vtxPassID')   


    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    wgt = events.eventWgt/sum_wgt
    wgt_frac = events.genWgt/sum_wgt
    
    if info["type"] == "signal":
        hasMatchVtx = ak.values_astype(ak.any(events.vtx.isMatched,axis=1),int)
        signalReconstructed = ak.values_astype(events.signalReconstructed,int)
        signalReconstructedPassID = ak.values_astype(events.GenEle.matchPassID & events.GenPos.matchPassID,int)
        signalReconstructedPassIDBasic = ak.values_astype(events.GenEle.matchPassIDBasic & events.GenPos.matchPassIDBasic,int)
        matchVtxPassID = ak.values_astype(events.matchedVtxPassID,int)
        matchVtxPassIDBasic = ak.values_astype(events.matchedVtxPassIDBasic,int)

        eMatch = ak.values_astype(events.GenEle.matched,int)
        eMatchPassID = ak.values_astype(events.GenEle.matchPassID,int)
        eMatchPassIDBasic = ak.values_astype(events.GenEle.matchPassIDBasic,int)

        pMatch = ak.values_astype(events.GenPos.matched,int)
        pMatchPassID = ak.values_astype(events.GenPos.matchPassID,int)
        pMatchPassIDBasic = ak.values_astype(events.GenPos.matchPassIDBasic,int)
        ### FILLING HISTOGRAMS ###
        
        #
        h.fill('signalReco_vs_vtxMatch_BaseID',sigReco=signalReconstructed,vtxMatch=hasMatchVtx,sigRecoPassID=signalReconstructedPassIDBasic,vtxPassID=matchVtxPassIDBasic,weight=wgt_frac)
        h.fill('signalReco_vs_vtxMatch_dRID',sigReco=signalReconstructed,vtxMatch=hasMatchVtx,sigRecoPassID=signalReconstructedPassID,vtxPassID=matchVtxPassID,weight=wgt_frac)
        h.fill('genEPMatched_BaseID',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassIDBasic,pID=pMatchPassIDBasic,weight=wgt_frac)
        h.fill('genEPMatched_dRID',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassID,pID=pMatchPassID,weight=wgt_frac)
        h.fill('genEPMatched_vtxMatched_BaseID',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassIDBasic,pID=pMatchPassIDBasic,vtxMatch=hasMatchVtx,vtxPassID=matchVtxPassIDBasic,weight=wgt_frac)
        h.fill('genEPMatched_vtxMatched_dRID',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassID,pID=pMatchPassID,vtxMatch=hasMatchVtx,vtxPassID=matchVtxPassID,weight=wgt_frac)

        h.fill('signalReco_vs_vtxMatch_BaseID_unwgt',sigReco=signalReconstructed,vtxMatch=hasMatchVtx,sigRecoPassID=signalReconstructedPassIDBasic,vtxPassID=matchVtxPassIDBasic,weight=1)
        h.fill('signalReco_vs_vtxMatch_dRID_unwgt',sigReco=signalReconstructed,vtxMatch=hasMatchVtx,sigRecoPassID=signalReconstructedPassID,vtxPassID=matchVtxPassID,weight=1)
        h.fill('genEPMatched_BaseID_unwgt',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassIDBasic,pID=pMatchPassIDBasic,weight=1)
        h.fill('genEPMatched_dRID_unwgt',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassID,pID=pMatchPassID,weight=1)
        h.fill('genEPMatched_vtxMatched_BaseID_unwgt',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassIDBasic,pID=pMatchPassIDBasic,vtxMatch=hasMatchVtx,vtxPassID=matchVtxPassIDBasic,weight=1)
        h.fill('genEPMatched_vtxMatched_dRID_unwgt',eMatch=eMatch,pMatch=pMatch,eID=eMatchPassID,pID=pMatchPassID,vtxMatch=hasMatchVtx,vtxPassID=matchVtxPassID,weight=1)
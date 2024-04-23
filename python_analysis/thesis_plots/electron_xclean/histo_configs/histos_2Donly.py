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
    h.make('gen_ele_pt','ele_pt')
    h.make('gen_ele_eta','eta')
    h.make('gen_ele_phi','phi')
    h.make('gen_ele_vxy1','vxy1')
    h.make('gen_ele_vxy10','vxy10')
    h.make('gen_ele_vxy100','vxy100')
    h.make("gen_ele_pt_vs_vxy1",'ele_pt','vxy1')
    h.make("gen_ele_pt_vs_vxy10",'ele_pt','vxy10')
    h.make("gen_ele_pt_vs_vxy100",'ele_pt','vxy100')
    h.make("gen_ele_pt_vs_dR",'ele_pt','dR')
    h.make("gen_ele_dR_vs_vxy1",'dR','vxy1')
    h.make("gen_ele_dR_vs_vxy10",'dR','vxy10')
    h.make("gen_ele_dR_vs_vxy100",'dR','vxy100')
    h.make("gen_ele_pt_vs_vxy1_dRBinned",'ele_pt','vxy1','dRCategories')
    h.make("gen_ele_pt_vs_vxy10_dRBinned",'ele_pt','vxy10','dRCategories')
    h.make("gen_ele_pt_vs_vxy100_dRBinned",'ele_pt','vxy100','dRCategories')
    h.make("gen_ele_pt_vs_dR_vxyBinned",'ele_pt','dR','vxyCategories')
    h.make("gen_ele_dR_vs_vxy1_ptBinned",'dR','vxy1','ptCategories')
    h.make("gen_ele_dR_vs_vxy10_ptBinned",'dR','vxy10','ptCategories')
    h.make("gen_ele_dR_vs_vxy100_ptBinned",'dR','vxy100','ptCategories')
    
    for etype in ['L','R','Both']:
        for id in [0,1,2,3]:
            # matched reco electrons, corresponding gen object variables, save ID result including dR(e,jets) cut
            h.make(f"match_ele_gen_pt_Type{etype}_ID{id}",'ele_pt')
            h.make(f'match_ele_gen_vxy1_Type{etype}_ID{id}','vxy1')
            h.make(f'match_ele_gen_vxy10_Type{etype}_ID{id}','vxy10')
            h.make(f'match_ele_gen_vxy100_Type{etype}_ID{id}','vxy100')
            h.make(f'match_ele_gen_pt_vs_vxy1_Type{etype}_ID{id}','ele_pt','vxy1')
            h.make(f'match_ele_gen_pt_vs_vxy10_Type{etype}_ID{id}','ele_pt','vxy10')
            h.make(f'match_ele_gen_pt_vs_vxy100_Type{etype}_ID{id}','ele_pt','vxy100')
            h.make(f"match_ele_gen_pt_vs_dR_Type{etype}_ID{id}",'ele_pt','dR')
            h.make(f"match_ele_gen_dR_vs_vxy1_Type{etype}_ID{id}",'dR','vxy1')
            h.make(f"match_ele_gen_dR_vs_vxy10_Type{etype}_ID{id}",'dR','vxy10')
            h.make(f"match_ele_gen_dR_vs_vxy100_Type{etype}_ID{id}",'dR','vxy100')
            h.make(f"match_ele_gen_eta_Type{etype}_ID{id}",'eta')
            h.make(f'match_ele_gen_phi_Type{etype}_ID{id}','phi')
            h.make(f'match_ele_gen_dR_Type{etype}_ID{id}','dR')
            h.make(f"match_ele_gen_pt_vs_vxy1_dRBinned_Type{etype}_ID{id}",'ele_pt','vxy1','dRCategories')
            h.make(f"match_ele_gen_pt_vs_vxy10_dRBinned_Type{etype}_ID{id}",'ele_pt','vxy10','dRCategories')
            h.make(f"match_ele_gen_pt_vs_vxy100_dRBinned_Type{etype}_ID{id}",'ele_pt','vxy100','dRCategories')
            h.make(f"match_ele_gen_pt_vs_dR_vxyBinned_Type{etype}_ID{id}",'ele_pt','dR','vxyCategories')
            h.make(f"match_ele_gen_dR_vs_vxy1_ptBinned_Type{etype}_ID{id}",'dR','vxy1','ptCategories')
            h.make(f"match_ele_gen_dR_vs_vxy10_ptBinned_Type{etype}_ID{id}",'dR','vxy10','ptCategories')
            h.make(f"match_ele_gen_dR_vs_vxy100_ptBinned_Type{etype}_ID{id}",'dR','vxy100','ptCategories')

            # matched reco electrons, corresponding gen object variables, save ID result without dR(e,jets) cut
            h.make(f"match_ele_gen_pt_IDbasic_Type{etype}_ID{id}",'ele_pt')
            h.make(f'match_ele_gen_vxy1_IDbasic_Type{etype}_ID{id}','vxy1')
            h.make(f'match_ele_gen_vxy10_IDbasic_Type{etype}_ID{id}','vxy10')
            h.make(f'match_ele_gen_vxy100_IDbasic_Type{etype}_ID{id}','vxy100')
            h.make(f'match_ele_gen_pt_vs_vxy1_IDbasic_Type{etype}_ID{id}','ele_pt','vxy1')
            h.make(f'match_ele_gen_pt_vs_vxy10_IDbasic_Type{etype}_ID{id}','ele_pt','vxy10')
            h.make(f'match_ele_gen_pt_vs_vxy100_IDbasic_Type{etype}_ID{id}','ele_pt','vxy100')
            h.make(f"match_ele_gen_pt_vs_dR_IDbasic_Type{etype}_ID{id}",'ele_pt','dR')
            h.make(f"match_ele_gen_dR_vs_vxy1_IDbasic_Type{etype}_ID{id}",'dR','vxy1')
            h.make(f"match_ele_gen_dR_vs_vxy10_IDbasic_Type{etype}_ID{id}",'dR','vxy10')
            h.make(f"match_ele_gen_dR_vs_vxy100_IDbasic_Type{etype}_ID{id}",'dR','vxy100')
            h.make(f"match_ele_gen_eta_IDbasic_Type{etype}_ID{id}",'eta')
            h.make(f'match_ele_gen_phi_IDbasic_Type{etype}_ID{id}','phi')
            h.make(f'match_ele_gen_dR_IDbasic_Type{etype}_ID{id}','dR')
            h.make(f"match_ele_gen_pt_vs_vxy1_dRBinned_IDbasic_Type{etype}_ID{id}",'ele_pt','vxy1','dRCategories')
            h.make(f"match_ele_gen_pt_vs_vxy10_dRBinned_IDbasic_Type{etype}_ID{id}",'ele_pt','vxy10','dRCategories')
            h.make(f"match_ele_gen_pt_vs_vxy100_dRBinned_IDbasic_Type{etype}_ID{id}",'ele_pt','vxy100','dRCategories')
            h.make(f"match_ele_gen_pt_vs_dR_vxyBinned_IDbasic_Type{etype}_ID{id}",'ele_pt','dR','vxyCategories')
            h.make(f"match_ele_gen_dR_vs_vxy1_ptBinned_IDbasic_Type{etype}_ID{id}",'dR','vxy1','ptCategories')
            h.make(f"match_ele_gen_dR_vs_vxy10_ptBinned_IDbasic_Type{etype}_ID{id}",'dR','vxy10','ptCategories')
            h.make(f"match_ele_gen_dR_vs_vxy100_ptBinned_IDbasic_Type{etype}_ID{id}",'dR','vxy100','ptCategories')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    wgt = events.eventWgt/sum_wgt
    
    if info["type"] == "signal":
        # defining stuff
        hasMatch_pf = (ak.count(events.Electron.pt,axis=1)>0) &\
                      (ak.count_nonzero(events.Electron.genMatched & ~events.Electron.hasLptMatch,axis=1)>0)

        hasMatch_both = (ak.count(events.Electron.pt,axis=1)>0) &\
                        (ak.count_nonzero(events.Electron.genMatched & events.Electron.hasLptMatch,axis=1)>0)

        hasMatch_lpt = (ak.count(events.LptElectron.pt,axis=1)>0) &\
                       (ak.count_nonzero(events.LptElectron.genMatched,axis=1)>0)

        match_pf = events[hasMatch_pf].Electron
        match_pf = match_pf[match_pf.genMatched & ~match_pf.hasLptMatch]
        genObj_pf = ak.where(match_pf.matchType==-1,events[hasMatch_pf].GenEle,events[hasMatch_pf].GenPos)
        match_pf = ak.flatten(match_pf)
        genObj_pf = ak.flatten(genObj_pf)
        match_pf_passID = ak.values_astype(match_pf.passID,int)
        match_pf_passIDBasic = ak.values_astype(match_pf.passIDBasic,int)

        match_both = events[hasMatch_both].Electron
        match_both = match_both[match_both.genMatched & match_both.hasLptMatch]
        genObj_both = ak.where(match_both.matchType==-1,events[hasMatch_both].GenEle,events[hasMatch_both].GenPos)
        lptObj_both = ak.flatten(events[hasMatch_both].LptElectron[match_both.lptMatchIdx])
        match_both = ak.flatten(match_both)
        genObj_both = ak.flatten(genObj_both)
        match_both_passID = ak.zeros_like(match_both.pt,dtype=np.int32)
        match_both_passID = ak.where(match_both.passID & lptObj_both.passID,1,match_both_passID)
        match_both_passID = ak.where(match_both.passID & ~lptObj_both.passID,2,match_both_passID)
        match_both_passID = ak.where(~match_both.passID & lptObj_both.passID,3,match_both_passID)
        match_both_passIDBasic = ak.zeros_like(match_both.pt,dtype=np.int32)
        match_both_passIDBasic = ak.where(match_both.passIDBasic & lptObj_both.passIDBasic,1,match_both_passIDBasic)
        match_both_passIDBasic = ak.where(match_both.passIDBasic & ~lptObj_both.passIDBasic,2,match_both_passIDBasic)
        match_both_passIDBasic = ak.where(~match_both.passIDBasic & lptObj_both.passIDBasic,3,match_both_passIDBasic)

        match_lpt = events[hasMatch_lpt].LptElectron
        match_lpt = match_lpt[match_lpt.genMatched]
        genObj_lpt = ak.where(match_lpt.matchType==-1,events[hasMatch_lpt].GenEle,events[hasMatch_lpt].GenPos)
        match_lpt = ak.flatten(match_lpt)
        genObj_lpt = ak.flatten(genObj_lpt)
        match_lpt_passID = ak.values_astype(match_lpt.passID,int)
        match_lpt_passIDBasic = ak.values_astype(match_lpt.passIDBasic,int)

        matches = {
            "R":match_pf,
            "L":match_lpt,
            "Both":match_both
        }
        genObjs = {
            "R":genObj_pf,
            "L":genObj_lpt,
            "Both":genObj_both
        }
        passIDs = {
            "R":match_pf_passID,
            "L":match_lpt_passID,
            "Both":match_both_passID
        }
        passIDsBasic = {
            "R":match_pf_passIDBasic,
            "L":match_lpt_passIDBasic,
            "Both":match_both_passIDBasic
        }

        del hasMatch_pf, hasMatch_lpt, hasMatch_both

        ### FILLING HISTOGRAMS ###
        h.fill("gen_ele_pt",pt=events.GenEle.pt,weight=1)
        h.fill("gen_ele_pt",pt=events.GenPos.pt,weight=1)
        h.fill("gen_ele_vxy1",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy1",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_vxy10",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy10",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_vxy100",vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_vxy100",vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy1",pt=events.GenEle.pt,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy1",pt=events.GenPos.pt,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy10",pt=events.GenEle.pt,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy10",pt=events.GenPos.pt,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy100",pt=events.GenEle.pt,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy100",pt=events.GenPos.pt,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_pt_vs_dR",pt=events.GenEle.pt,dR=events.GenEle.dr,weight=1)
        h.fill("gen_ele_pt_vs_dR",pt=events.GenPos.pt,dR=events.GenPos.dr,weight=1)
        h.fill("gen_ele_dR_vs_vxy1",dR=events.GenEle.dr,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_dR_vs_vxy1",dR=events.GenPos.dr,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_dR_vs_vxy10",dR=events.GenEle.dr,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_dR_vs_vxy10",dR=events.GenPos.dr,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_dR_vs_vxy100",dR=events.GenEle.dr,vxy=events.GenEle.vxy,weight=1)
        h.fill("gen_ele_dR_vs_vxy100",dR=events.GenPos.dr,vxy=events.GenPos.vxy,weight=1)
        h.fill("gen_ele_pt_vs_vxy1_dRBinned",pt=events.GenEle.pt,vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_vxy1_dRBinned",pt=events.GenPos.pt,vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_vxy10_dRBinned",pt=events.GenEle.pt,vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_vxy10_dRBinned",pt=events.GenPos.pt,vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_vxy100_dRBinned",pt=events.GenEle.pt,vxy=events.GenEle.vxy,dRCat=events.GenEle.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_vxy100_dRBinned",pt=events.GenPos.pt,vxy=events.GenPos.vxy,dRCat=events.GenPos.dRbin,weight=1)
        h.fill("gen_ele_pt_vs_dR_vxyBinned",pt=events.GenEle.pt,dR=events.GenEle.dr,vxyCat=events.GenEle.vxyBin,weight=1)
        h.fill("gen_ele_pt_vs_dR_vxyBinned",pt=events.GenPos.pt,dR=events.GenPos.dr,vxyCat=events.GenPos.vxyBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy1_ptBinned",dR=events.GenEle.dr,vxy=events.GenEle.vxy,ptCat=events.GenEle.ptBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy1_ptBinned",dR=events.GenPos.dr,vxy=events.GenPos.vxy,ptCat=events.GenPos.ptBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy10_ptBinned",dR=events.GenEle.dr,vxy=events.GenEle.vxy,ptCat=events.GenEle.ptBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy10_ptBinned",dR=events.GenPos.dr,vxy=events.GenPos.vxy,ptCat=events.GenPos.ptBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy100_ptBinned",dR=events.GenEle.dr,vxy=events.GenEle.vxy,ptCat=events.GenEle.ptBin,weight=1)
        h.fill("gen_ele_dR_vs_vxy100_ptBinned",dR=events.GenPos.dr,vxy=events.GenPos.vxy,ptCat=events.GenPos.ptBin,weight=1)
        h.fill("gen_ele_eta",eta=events.GenEle.eta,weight=1)
        h.fill("gen_ele_eta",eta=events.GenPos.eta,weight=1)
        h.fill("gen_ele_phi",phi=events.GenEle.phi,weight=1)
        h.fill("gen_ele_phi",phi=events.GenPos.phi,weight=1)
        
        #
        for etype in ["L","R","Both"]:
            for id in [0,1,2,3]:
                targ_id = genObjs[etype][passIDs[etype]==id]
                targ_idbasic = genObjs[etype][passIDsBasic[etype]==id]

                h.fill(f"match_ele_gen_pt_Type{etype}_ID{id}",pt=targ_id.pt,weight=1)
                h.fill(f"match_ele_gen_vxy1_Type{etype}_ID{id}",vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_vxy10_Type{etype}_ID{id}",vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_vxy100_Type{etype}_ID{id}",vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy1_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy10_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy100_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_dR_Type{etype}_ID{id}",pt=targ_id.pt,dR=targ_id.dr,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy1_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy10_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy100_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,weight=1)
                h.fill(f"match_ele_gen_eta_Type{etype}_ID{id}",eta=targ_id.eta,weight=1)
                h.fill(f"match_ele_gen_phi_Type{etype}_ID{id}",phi=targ_id.phi,weight=1)
                h.fill(f"match_ele_gen_dR_Type{etype}_ID{id}",dR=targ_id.dr,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy1_dRBinned_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,dRCat=targ_id.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy10_dRBinned_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,dRCat=targ_id.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy100_dRBinned_Type{etype}_ID{id}",pt=targ_id.pt,vxy=targ_id.vxy,dRCat=targ_id.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_dR_vxyBinned_Type{etype}_ID{id}",pt=targ_id.pt,dR=targ_id.dr,vxyCat=targ_id.vxyBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy1_ptBinned_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,ptCat=targ_id.ptBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy10_ptBinned_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,ptCat=targ_id.ptBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy100_ptBinned_Type{etype}_ID{id}",dR=targ_id.dr,vxy=targ_id.vxy,ptCat=targ_id.ptBin,weight=1)

                #
                h.fill(f"match_ele_gen_pt_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,weight=1)
                h.fill(f"match_ele_gen_vxy1_IDbasic_Type{etype}_ID{id}",
                        vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_vxy10_IDbasic_Type{etype}_ID{id}",
                        vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_vxy100_IDbasic_Type{etype}_ID{id}",
                        vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy1_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy10_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy100_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_pt_vs_dR_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,dR=targ_idbasic.dr,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy1_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy10_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy100_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,weight=1)
                h.fill(f"match_ele_gen_eta_IDbasic_Type{etype}_ID{id}",
                        eta=targ_idbasic.eta,weight=1)
                h.fill(f"match_ele_gen_phi_IDbasic_Type{etype}_ID{id}",
                        phi=targ_idbasic.phi,weight=1)
                h.fill(f"match_ele_gen_dR_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy1_dRBinned_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,dRCat=targ_idbasic.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy10_dRBinned_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,dRCat=targ_idbasic.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_vxy100_dRBinned_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,vxy=targ_idbasic.vxy,dRCat=targ_idbasic.dRbin,weight=1)
                h.fill(f"match_ele_gen_pt_vs_dR_vxyBinned_IDbasic_Type{etype}_ID{id}",
                        pt=targ_idbasic.pt,dR=targ_idbasic.dr,vxyCat=targ_idbasic.vxyBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy1_ptBinned_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,ptCat=targ_idbasic.ptBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy10_ptBinned_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,ptCat=targ_idbasic.ptBin,weight=1)
                h.fill(f"match_ele_gen_dR_vs_vxy100_ptBinned_IDbasic_Type{etype}_ID{id}",
                        dR=targ_idbasic.dr,vxy=targ_idbasic.vxy,ptCat=targ_idbasic.ptBin,weight=1)
                
                del targ_id, targ_idbasic
        
        del match_both,match_lpt,match_pf
        del genObj_both,genObj_lpt,genObj_pf
        del match_pf_passID, match_lpt_passID, match_both_passID
        del match_pf_passIDBasic,match_lpt_passIDBasic,match_both_passIDBasic
        del matches,genObjs,passIDs,passIDsBasic
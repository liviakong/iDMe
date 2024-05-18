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
        self.abs_dphi = self.parse_axis(('abs_dphi',100,0,3.2))
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
        self.nJets = self.parse_axis(('nJets',10,0,10))
        self.mass_low = self.parse_axis(('mass_low',500,0,5))
        self.mindxy_low = self.parse_axis(('mindxy_low',100,0,0.1))
        self.sign_etaProd = self.parse_axis(('sign_etaProd',[-1,1]))
        self.cosTheta = self.parse_axis(('cosTheta',100,-1,1))
        self.LxyCosTheta = self.parse_axis(('LxyCosTheta',100,-50,50))
        self.LxyCosThetaZoom = self.parse_axis(('LxyCosThetaZoom',100,-5,5))
        self.LxyCosThetaZoomZoom = self.parse_axis(('LxyCosThetaZoomZoom',100,-1,1))
        self.jetMETratio = self.parse_axis(('jetMETratio',100,0,2))
        self.vtxPurity = self.parse_axis(('purity',['total','n_matched','n_eeReco','n_vtxReco']))
        self.chi2Rank = self.parse_axis(('chi2Rank',[0,1,2,3,4,5,6,7,8,9,10]))

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
    h.make("gen_vtx_METdPhi",'abs_dphi')
    h.make("gen_jetMETdPhi",'abs_dphi')
    h.make("gen_vtx_mass",'vtx_mass')
    h.make('gen_vtx_pt','vtx_pt')
    h.make('gen_vtx_eta','eta')
    h.make('gen_vtx_phi','phi')

    h.make('gen_ele_pt','ele_pt')
    h.make('gen_ele_eta','eta')
    h.make('gen_ele_phi','phi')
    h.make('gen_ele_vxy1','vxy1')
    h.make('gen_ele_vxy10','vxy10')
    h.make('gen_ele_vxy100','vxy100')
    h.make("gen_ele_r3", 'r3')
    h.make("gen_ele_r3_PVcorr", 'r3')

    # reco --> vertexing efficiency
    h.make("signalReco_vs_vtxMatch",'sigReco','vtxMatch')
    h.make("signalReco_vs_vtxMatch_unwgt",'sigReco','vtxMatch')
    
    # matched reco electron, reco variables
    h.make("match_ele_pt",'match_type','ele_pt')
    h.make("match_ele_eta",'match_type','eta')
    h.make("match_ele_phi",'match_type','phi')
    h.make("match_ele_dxy",'match_type','dxy')
    h.make("match_ele_dxySignif",'match_type','dxy_signif')
    h.make("match_ele_trkChi2",'match_type','trk_chi2')
    h.make('match_ele_trkProb','match_type','trkProb')
    h.make("match_ele_trkRelIso",'match_type','iso')
    h.make("match_ele_calRelIso",'match_type','iso')
    h.make("match_ele_PFRelIso",'match_type','iso')
    h.make("match_ele_miniRelIso",'match_type','iso')
    h.make("match_ele_mindRJets",'match_type','dR')
    h.make("match_ele_isPF",'match_type','isPF')
    h.make("match_ele_numTrkHits",'match_type','numHits')
    h.make("match_ele_numPixHits",'match_type','numHits')
    h.make("match_ele_numStripHits",'match_type','numHits')
    h.make("match_ele_mindPhiJets",'match_type','abs_dphi')
    h.make("match_ele_sigmaIetaIeta",'match_type','sieie')
    h.make("match_ele_absdEtaSeed",'match_type','detaseed')
    h.make("match_ele_absdPhiIn",'match_type','abs_dphi')
    h.make("match_ele_HOverE",'match_type','HoverE')
    h.make("match_ele_1E1p",'match_type','EinvMinusPinv')
    h.make("match_ele_expMissing",'match_type','numMissingHits')
    h.make("match_ele_passConvVeto",'match_type','passConvVeto')
    h.make("match_ele_IDscore",'match_type','IDScore')

    # fake reco electrons (not matched to a gen object)
    h.make("fake_ele_pt",'ele_type','ele_pt')
    h.make("fake_ele_eta",'ele_type','eta')
    h.make("fake_ele_phi",'ele_type','phi')
    h.make("fake_ele_dxy",'ele_type','dxy')
    h.make("fake_ele_dxySignif",'ele_type','dxy_signif')
    h.make("fake_ele_trkChi2",'ele_type','trk_chi2')
    h.make('fake_ele_trkProb','ele_type','trkProb')
    h.make("fake_ele_trkRelIso",'ele_type','iso')
    h.make("fake_ele_calRelIso",'ele_type','iso')
    h.make("fake_ele_PFRelIso",'ele_type','iso')
    h.make("fake_ele_miniRelIso",'ele_type','iso')
    h.make("fake_ele_mindRJets",'ele_type','dR')
    h.make("fake_ele_isPF",'ele_type','isPF')
    h.make("fake_ele_numTrkHits",'ele_type','numHits')
    h.make("fake_ele_numPixHits",'ele_type','numHits')
    h.make("fake_ele_numStripHits",'ele_type','numHits')
    h.make("fake_ele_mindPhiJets",'ele_type','abs_dphi')
    h.make("fake_ele_sigmaIetaIeta",'ele_type','sieie')
    h.make("fake_ele_absdEtaSeed",'ele_type','detaseed')
    h.make("fake_ele_absdPhiIn",'ele_type','abs_dphi')
    h.make("fake_ele_HOverE",'ele_type','HoverE')
    h.make("fake_ele_1E1p",'ele_type','EinvMinusPinv')
    h.make("fake_ele_expMissing",'ele_type','numMissingHits')
    h.make("fake_ele_passConvVeto",'ele_type','passConvVeto')
    h.make("fake_ele_IDscore",'ele_type','IDScore')
    
    # matched reco vertex
    h.make("match_vtx_dR",'vtx_type','dR')
    h.make("match_vtx_mindxy",'vtx_type','dxy')
    h.make("match_vtx_vxy1",'vtx_type','vxy1')
    h.make("match_vtx_vxy10",'vtx_type','vxy10')
    h.make("match_vtx_vxy100",'vtx_type','vxy100')
    h.make("match_vtx_leadpT",'vtx_type','ele_pt')
    h.make("match_vtx_METdPhi",'vtx_type','abs_dphi')
    h.make("match_vtx_mindRj",'vtx_type','dR')
    h.make("match_vtx_chi2",'vtx_type','vtx_chi2')
    h.make('match_vtx_mass','vtx_type','vtx_mass')
    h.make('match_vtx_mindPhiJ','vtx_type','abs_dphi')
    h.make('match_vtx_sign','vtx_type','vtx_sign')
    h.make('match_vtx_pt','vtx_type','vtx_pt')
    h.make('match_vtx_eta','vtx_type','eta')
    h.make('match_vtx_phi','vtx_type','phi')
    h.make("match_vtx_type",'vtx_type')
    h.make("match_vtx_minEleDrJ",'vtx_type','dR')
    h.make("match_vtx_minEleDPhiJ",'vtx_type','abs_dphi')
    h.make("match_vtx_mass_low",'vtx_type','mass_low')
    h.make("match_vtx_mindxy_low",'vtx_type','mindxy_low')
    h.make("match_vtx_sign_etaProd",'vtx_type','sign_etaProd')
    h.make("match_vtx_CosThetaColl",'vtx_type','cosTheta')
    h.make("match_vtx_LxyCosThetaColl",'vtx_type','LxyCosTheta')
    h.make("match_vtx_LxyCosThetaCollZoom",'vtx_type','LxyCosThetaZoom')
    h.make("match_vtx_LxyCosThetaCollZoomZoom",'vtx_type','LxyCosThetaZoomZoom')
    h.make("match_vtx_eleDphi",'vtx_type','abs_dphi')
    h.make("match_vtx_chi2Rank",'vtx_type','chi2Rank')

    # matched reco vertex, corresponding gen features
    h.make("match_vtx_gen_dR",'vtx_type','dR')
    h.make("match_vtx_gen_vxy1",'vtx_type','vxy1')
    h.make("match_vtx_gen_vxy10",'vtx_type','vxy10')
    h.make("match_vtx_gen_vxy100",'vtx_type','vxy100')
    h.make("match_vtx_gen_METdPhi",'vtx_type','abs_dphi')
    h.make("match_vtx_gen_mass",'vtx_type','vtx_mass')
    h.make('match_vtx_gen_pt','vtx_type','vtx_pt')
    h.make('match_vtx_gen_eta','vtx_type','eta')
    h.make('match_vtx_gen_phi','vtx_type','phi')

    # fake reco vertices
    h.make("fake_vtx_dR",'vtx_type','dR')
    h.make("fake_vtx_mindxy",'vtx_type','dxy')
    h.make("fake_vtx_vxy1",'vtx_type','vxy1')
    h.make("fake_vtx_vxy10",'vtx_type','vxy10')
    h.make("fake_vtx_vxy100",'vtx_type','vxy100')
    h.make("fake_vtx_leadpT",'vtx_type','ele_pt')
    h.make("fake_vtx_METdPhi",'vtx_type','abs_dphi')
    h.make("fake_vtx_mindRj",'vtx_type','dR')
    h.make("fake_vtx_chi2",'vtx_type','vtx_chi2')
    h.make('fake_vtx_mass','vtx_type','vtx_mass')
    h.make('fake_vtx_mindPhiJ','vtx_type','abs_dphi')
    h.make('fake_vtx_sign','vtx_type','vtx_sign')
    h.make('fake_vtx_pt','vtx_type','vtx_pt')
    h.make('fake_vtx_eta','vtx_type','eta')
    h.make('fake_vtx_phi','vtx_type','phi')
    h.make("fake_vtx_type",'vtx_type')
    h.make("fake_vtx_minEleDrJ",'vtx_type','dR')
    h.make("fake_vtx_minEleDPhiJ",'vtx_type','abs_dphi')
    h.make("fake_vtx_mass_low",'vtx_type','mass_low')
    h.make("fake_vtx_mindxy_low",'vtx_type','mindxy_low')
    h.make("fake_vtx_sign_etaProd",'vtx_type','sign_etaProd')
    h.make("fake_vtx_CosThetaColl",'vtx_type','cosTheta')
    h.make("fake_vtx_LxyCosThetaColl",'vtx_type','LxyCosTheta')
    h.make("fake_vtx_LxyCosThetaCollZoom",'vtx_type','LxyCosThetaZoom')
    h.make("fake_vtx_LxyCosThetaCollZoomZoom",'vtx_type','LxyCosThetaZoomZoom')
    h.make("fake_vtx_eleDphi",'vtx_type','abs_dphi')
    
    # misc other quantities
    h.make("PFMET",'met')
    h.make("PFMET_vs_genMET",'met',('genmet',100,50,300))
    h.make("jetMETdPhi",'abs_dphi')
    h.make("minJetMETdPhi",'abs_dphi')
    h.make("nJets",'nJets')
    h.make("jetMETratio",'jetMETratio')
    h.make('sel_vtx_isMatched','isMatched')
    h.make('sel_vtx_purity','vtxPurity')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    wgt = events.eventWgt/sum_wgt
    
    if info["type"] == "signal":
        # defining stuff
        hasMatch_pf = (ak.count(events.Electron.pt,axis=1)>0) &\
                      (ak.count_nonzero(events.Electron.genMatched & ~events.Electron.hasLptMatch & events.Electron.passID,axis=1)>0)

        hasMatch_both = (ak.count(events.Electron.pt,axis=1)>0) &\
                        (ak.count_nonzero(events.Electron.genMatched & events.Electron.hasLptMatch,axis=1)>0)

        hasMatch_lpt = (ak.count(events.LptElectron.pt,axis=1)>0) &\
                       (ak.count_nonzero(events.LptElectron.genMatched & events.LptElectron.passID,axis=1)>0)

        match_pf = events[hasMatch_pf].Electron
        match_pf = match_pf[match_pf.genMatched & ~match_pf.hasLptMatch & match_pf.passID]
        genObj_pf = ak.where(match_pf.matchType==-1,events[hasMatch_pf].GenEle,events[hasMatch_pf].GenPos)
        match_pf = ak.flatten(match_pf)
        genObj_pf = ak.flatten(genObj_pf)

        match_both = events[hasMatch_both].Electron
        match_both = match_both[match_both.genMatched & match_both.hasLptMatch]
        lptObj_both = events[hasMatch_both].LptElectron[match_both.lptMatchIdx]
        cut_passIDAny_both = ak.any(match_both.passID | lptObj_both.passID,axis=1)
        match_both = match_both[cut_passIDAny_both]
        lptObj_both = lptObj_both[cut_passIDAny_both]
        cut_selPassID_both = match_both.passID | lptObj_both.passID
        match_both = match_both[cut_selPassID_both]
        lptObj_both = lptObj_both[cut_selPassID_both]
        genObj_both = ak.where(match_both.matchType==-1,events[hasMatch_both][cut_passIDAny_both].GenEle,events[hasMatch_both][cut_passIDAny_both].GenPos)
        match_both = ak.flatten(match_both)
        lptObj_both = ak.flatten(lptObj_both)
        genObj_both = ak.flatten(genObj_both)

        match_lpt = events[hasMatch_lpt].LptElectron
        match_lpt = match_lpt[match_lpt.genMatched & match_lpt.passID]
        genObj_lpt = ak.where(match_lpt.matchType==-1,events[hasMatch_lpt].GenEle,events[hasMatch_lpt].GenPos)
        match_lpt = ak.flatten(match_lpt)
        genObj_lpt = ak.flatten(genObj_lpt)

        fake_pf = ak.flatten(events.Electron[~events.Electron.genMatched & ~events.Electron.hasLptMatch & events.Electron.passID])
        
        hasfake_both = ak.any(~events.Electron.genMatched & events.Electron.hasLptMatch,axis=1)
        fake_both = events[hasfake_both].Electron
        fake_both = fake_both[~fake_both.genMatched & fake_both.hasLptMatch]
        lptObj_fake_both = events[hasfake_both].LptElectron[fake_both.lptMatchIdx]
        fake_both_passID_cut = ak.any(fake_both.passID | lptObj_fake_both.passID,axis=1)
        fake_both = fake_both[fake_both_passID_cut]
        lptObj_fake_both = lptObj_fake_both[fake_both_passID_cut]
        fake_both_passID_ele = fake_both.passID | lptObj_fake_both.passID
        fake_both = ak.flatten(fake_both[fake_both_passID_ele])
        
        fake_lpt = ak.flatten(events.LptElectron[~events.LptElectron.genMatched & ~events.LptElectron.gedIsMatched & events.LptElectron.passID])

        match_vtx_cut = events.good_vtx.isMatched
        fake_vtx_cut = ~events.good_vtx.isMatched
        match_vtx = ak.flatten(events.good_vtx[match_vtx_cut])
        match_vtx_genObj = events.genEE[ak.any(match_vtx_cut,axis=1)]
        match_vtx_wgt = wgt[ak.any(match_vtx_cut,axis=1)]
        fake_vtx = ak.flatten(events.good_vtx[fake_vtx_cut])

        hasMatch_vtx = ak.count_nonzero(events.good_vtx.isMatched,axis=1) > 0
        match_vtx_chi2Sort = ak.argsort(events[hasMatch_vtx].good_vtx.reduced_chi2,axis=1)
        match_vtx_chi2Rank = ak.flatten(match_vtx_chi2Sort[events[hasMatch_vtx].good_vtx.isMatched])

        hasMatchVtx = ak.values_astype(ak.any(events.vtx.isMatched,axis=1),int)
        signalReconstructed = ak.values_astype(events.signalReconstructed,int)

        gen_r3 = np.sqrt(events.GenEle.vx**2 + events.GenEle.vy**2 + events.GenEle.vz**2)
        gen_r3PV = np.sqrt((events.GenEle.vx - events.PV.x)**2 + (events.GenEle.vy - events.PV.y)**2 + (events.GenEle.vz - events.PV.z)**2)

        match_vtx_cosTheta = (match_vtx.px/match_vtx.pt)*(match_vtx.vx/match_vtx.vxy) + (match_vtx.py/match_vtx.pt)*(match_vtx.vy/match_vtx.vxy)
        fake_vtx_cosTheta = (fake_vtx.px/fake_vtx.pt)*(fake_vtx.vx/fake_vtx.vxy) + (fake_vtx.py/fake_vtx.pt)*(fake_vtx.vy/fake_vtx.vxy)

        ### FILLING HISTOGRAMS ###

        h.fill("gen_met",met=events.GenMET.pt,weight=wgt)
        h.fill("gen_dR",dR=events.genEE.dr,weight=wgt)
        h.fill("gen_vxy1",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_vxy10",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_vxy100",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_leadpT",pt=np.maximum(events.GenEle.pt,events.GenPos.pt),weight=wgt)
        h.fill("gen_vtx_METdPhi",abs_dphi=np.abs(events.genEE.METdPhi),weight=wgt)
        h.fill("gen_jetMETdPhi",abs_dphi=np.abs(events.GenJetMETdPhi),weight=wgt)
        h.fill("gen_vtx_mass",mass=events.genEE.mass,weight=wgt)
        h.fill('gen_vtx_pt',pt=events.genEE.pt,weight=wgt)
        h.fill('gen_vtx_eta',eta=events.genEE.eta,weight=wgt)
        h.fill('gen_vtx_phi',phi=events.genEE.phi,weight=wgt)
        
        h.fill("gen_ele_pt",pt=events.GenEle.pt,weight=wgt)
        h.fill("gen_ele_pt",pt=events.GenPos.pt,weight=wgt)
        h.fill("gen_ele_vxy1",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_ele_vxy1",vxy=events.GenPos.vxy,weight=wgt)
        h.fill("gen_ele_vxy10",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_ele_vxy10",vxy=events.GenPos.vxy,weight=wgt)
        h.fill("gen_ele_vxy100",vxy=events.GenEle.vxy,weight=wgt)
        h.fill("gen_ele_vxy100",vxy=events.GenPos.vxy,weight=wgt)
        h.fill("gen_ele_eta",eta=events.GenEle.eta,weight=wgt)
        h.fill("gen_ele_eta",eta=events.GenPos.eta,weight=wgt)
        h.fill("gen_ele_phi",phi=events.GenEle.phi,weight=wgt)
        h.fill("gen_ele_phi",phi=events.GenPos.phi,weight=wgt)
        h.fill("gen_ele_r3",r3=gen_r3,weight=wgt)
        h.fill("gen_ele_r3_PVcorr",r3=gen_r3PV,weight=wgt)
        
        #
        h.fill('signalReco_vs_vtxMatch',reco=signalReconstructed,match=hasMatchVtx,weight=wgt)
        h.fill('signalReco_vs_vtxMatch_unwgt',reco=signalReconstructed,match=hasMatchVtx,weight=1)  
        
        #
        h.fill("match_ele_pt",match_type='R',pt=match_pf.pt,weight=1)
        h.fill("match_ele_eta",match_type='R',eta=match_pf.eta,weight=1)
        h.fill("match_ele_phi",match_type='R',phi=match_pf.phi,weight=1)
        h.fill("match_ele_dxy",match_type='R',dxy=match_pf.dxy,weight=1)
        h.fill("match_ele_dxySignif",match_type='R',signif=np.abs(match_pf.dxy)/match_pf.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",match_type='R',chi2=match_pf.trkChi2,weight=1)
        h.fill("match_ele_trkProb",match_type='R',prob=match_pf.trkProb)
        h.fill("match_ele_trkRelIso",match_type='R',iso=match_pf.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",match_type='R',iso=match_pf.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",match_type='R',iso=match_pf.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",match_type='R',iso=match_pf.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",match_type='R',dR=match_pf.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",match_type='R',abs_dphi=np.abs(match_pf.mindPhiJ),weight=1)
        h.fill('match_ele_isPF',match_type='R',isPF=ak.values_astype(match_pf.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',match_type='R',numHits=match_pf.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',match_type='R',numHits=match_pf.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',match_type='R',numHits=match_pf.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",match_type='R',sieie=match_pf.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",match_type='R',detaSeed=match_pf.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",match_type='R',abs_dphi=np.abs(match_pf.absdPhiIn),weight=1)
        h.fill("match_ele_HOverE",match_type='R',hoe=match_pf.HoverE,weight=1)
        h.fill("match_ele_1E1p",match_type='R',emp=match_pf.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",match_type='R',missing=match_pf.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",match_type='R',passVeto=match_pf.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",match_type='R',id=match_pf.IDscore,weight=1)

        h.fill("match_ele_pt",match_type='Both',pt=match_both.pt,weight=1)
        h.fill("match_ele_eta",match_type='Both',eta=match_both.eta,weight=1)
        h.fill("match_ele_phi",match_type='Both',phi=match_both.phi,weight=1)
        h.fill("match_ele_dxy",match_type='Both',dxy=match_both.dxy,weight=1)
        h.fill("match_ele_dxySignif",match_type='Both',signif=np.abs(match_both.dxy)/match_both.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",match_type='Both',chi2=match_both.trkChi2,weight=1)
        h.fill("match_ele_trkProb",match_type='Both',prob=match_both.trkProb)
        h.fill("match_ele_trkRelIso",match_type='Both',iso=match_both.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",match_type='Both',iso=match_both.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",match_type='Both',iso=match_both.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",match_type='Both',iso=match_both.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",match_type='Both',dR=match_both.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",match_type='Both',abs_dphi=np.abs(match_both.mindPhiJ),weight=1)
        h.fill('match_ele_isPF',match_type='Both',isPF=ak.values_astype(match_both.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',match_type='Both',numHits=match_both.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',match_type='Both',numHits=match_both.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',match_type='Both',numHits=match_both.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",match_type='Both',sieie=match_both.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",match_type='Both',detaSeed=match_both.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",match_type='Both',abs_dphi=np.abs(match_both.absdPhiIn),weight=1)
        h.fill("match_ele_HOverE",match_type='Both',hoe=match_both.HoverE,weight=1)
        h.fill("match_ele_1E1p",match_type='Both',emp=match_both.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",match_type='Both',missing=match_both.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",match_type='Both',passVeto=match_both.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",match_type='Both',id=match_both.IDscore,weight=1)

        h.fill("match_ele_pt",match_type='L',pt=match_lpt.pt,weight=1)
        h.fill("match_ele_eta",match_type='L',eta=match_lpt.eta,weight=1)
        h.fill("match_ele_phi",match_type='L',phi=match_lpt.phi,weight=1)
        h.fill("match_ele_dxy",match_type='L',dxy=match_lpt.dxy,weight=1)
        h.fill("match_ele_dxySignif",match_type='L',signif=np.abs(match_lpt.dxy)/match_lpt.dxyErr,weight=1)
        h.fill("match_ele_trkChi2",match_type='L',chi2=match_lpt.trkChi2,weight=1)
        h.fill("match_ele_trkProb",match_type='L',prob=match_lpt.trkProb)
        h.fill("match_ele_trkRelIso",match_type='L',iso=match_lpt.trkRelIso,weight=1)
        h.fill("match_ele_calRelIso",match_type='L',iso=match_lpt.calRelIso,weight=1)
        h.fill("match_ele_PFRelIso",match_type='L',iso=match_lpt.PFRelIso,weight=1)
        h.fill("match_ele_miniRelIso",match_type='L',iso=match_lpt.miniRelIso,weight=1)
        h.fill("match_ele_mindRJets",match_type='L',dR=match_lpt.mindRj,weight=1)
        h.fill("match_ele_mindPhiJets",match_type='L',abs_dphi=np.abs(match_lpt.mindPhiJ),weight=1)
        h.fill('match_ele_isPF',match_type='L',isPF=ak.values_astype(match_lpt.isPF,int),weight=1)
        h.fill('match_ele_numTrkHits',match_type='L',numHits=match_lpt.numTrackerHits,weight=1)
        h.fill('match_ele_numPixHits',match_type='L',numHits=match_lpt.numPixHits,weight=1)
        h.fill('match_ele_numStripHits',match_type='L',numHits=match_lpt.numStripHits,weight=1)
        h.fill("match_ele_sigmaIetaIeta",match_type='L',sieie=match_lpt.full55sigmaIetaIeta,weight=1)
        h.fill("match_ele_absdEtaSeed",match_type='L',detaSeed=match_lpt.absdEtaSeed,weight=1)
        h.fill("match_ele_absdPhiIn",match_type='L',abs_dphi=np.abs(match_lpt.absdPhiIn),weight=1)
        h.fill("match_ele_HOverE",match_type='L',hoe=match_lpt.HoverE,weight=1)
        h.fill("match_ele_1E1p",match_type='L',emp=match_lpt.abs1overEm1overP,weight=1)
        h.fill("match_ele_expMissing",match_type='L',missing=match_lpt.expMissingInnerHits,weight=1)
        h.fill("match_ele_passConvVeto",match_type='L',passVeto=match_lpt.conversionVeto,weight=1)
        h.fill("match_ele_IDscore",match_type='L',id=match_lpt.IDscore,weight=1)

        #
        h.fill("fake_ele_pt",ele_type='R',pt=fake_pf.pt,weight=1)
        h.fill("fake_ele_eta",ele_type='R',eta=fake_pf.eta,weight=1)
        h.fill("fake_ele_phi",ele_type='R',phi=fake_pf.phi,weight=1)
        h.fill("fake_ele_dxy",ele_type='R',dxy=fake_pf.dxy,weight=1)
        h.fill("fake_ele_dxySignif",ele_type='R',signif=np.abs(fake_pf.dxy)/fake_pf.dxyErr,weight=1)
        h.fill("fake_ele_trkChi2",ele_type='R',chi2=fake_pf.trkChi2,weight=1)
        h.fill("fake_ele_trkProb",ele_type='R',prob=fake_pf.trkProb)
        h.fill("fake_ele_trkRelIso",ele_type='R',iso=fake_pf.trkRelIso,weight=1)
        h.fill("fake_ele_calRelIso",ele_type='R',iso=fake_pf.calRelIso,weight=1)
        h.fill("fake_ele_PFRelIso",ele_type='R',iso=fake_pf.PFRelIso,weight=1)
        h.fill("fake_ele_miniRelIso",ele_type='R',iso=fake_pf.miniRelIso,weight=1)
        h.fill("fake_ele_mindRJets",ele_type='R',dR=fake_pf.mindRj,weight=1)
        h.fill("fake_ele_mindPhiJets",ele_type='R',abs_dphi=np.abs(fake_pf.mindPhiJ),weight=1)
        h.fill('fake_ele_isPF',ele_type='R',isPF=ak.values_astype(fake_pf.isPF,int),weight=1)
        h.fill('fake_ele_numTrkHits',ele_type='R',numHits=fake_pf.numTrackerHits,weight=1)
        h.fill('fake_ele_numPixHits',ele_type='R',numHits=fake_pf.numPixHits,weight=1)
        h.fill('fake_ele_numStripHits',ele_type='R',numHits=fake_pf.numStripHits,weight=1)
        h.fill("fake_ele_sigmaIetaIeta",ele_type='R',sieie=fake_pf.full55sigmaIetaIeta,weight=1)
        h.fill("fake_ele_absdEtaSeed",ele_type='R',detaSeed=fake_pf.absdEtaSeed,weight=1)
        h.fill("fake_ele_absdPhiIn",ele_type='R',abs_dphi=np.abs(fake_pf.absdPhiIn),weight=1)
        h.fill("fake_ele_HOverE",ele_type='R',hoe=fake_pf.HoverE,weight=1)
        h.fill("fake_ele_1E1p",ele_type='R',emp=fake_pf.abs1overEm1overP,weight=1)
        h.fill("fake_ele_expMissing",ele_type='R',missing=fake_pf.expMissingInnerHits,weight=1)
        h.fill("fake_ele_passConvVeto",ele_type='R',passVeto=fake_pf.conversionVeto,weight=1)
        h.fill("fake_ele_IDscore",ele_type='R',id=fake_pf.IDscore,weight=1)

        h.fill("fake_ele_pt",ele_type='Both',pt=fake_both.pt,weight=1)
        h.fill("fake_ele_eta",ele_type='Both',eta=fake_both.eta,weight=1)
        h.fill("fake_ele_phi",ele_type='Both',phi=fake_both.phi,weight=1)
        h.fill("fake_ele_dxy",ele_type='Both',dxy=fake_both.dxy,weight=1)
        h.fill("fake_ele_dxySignif",ele_type='Both',signif=np.abs(fake_both.dxy)/fake_both.dxyErr,weight=1)
        h.fill("fake_ele_trkChi2",ele_type='Both',chi2=fake_both.trkChi2,weight=1)
        h.fill("fake_ele_trkProb",ele_type='Both',prob=fake_both.trkProb)
        h.fill("fake_ele_trkRelIso",ele_type='Both',iso=fake_both.trkRelIso,weight=1)
        h.fill("fake_ele_calRelIso",ele_type='Both',iso=fake_both.calRelIso,weight=1)
        h.fill("fake_ele_PFRelIso",ele_type='Both',iso=fake_both.PFRelIso,weight=1)
        h.fill("fake_ele_miniRelIso",ele_type='Both',iso=fake_both.miniRelIso,weight=1)
        h.fill("fake_ele_mindRJets",ele_type='Both',dR=fake_both.mindRj,weight=1)
        h.fill("fake_ele_mindPhiJets",ele_type='Both',abs_dphi=np.abs(fake_both.mindPhiJ),weight=1)
        h.fill('fake_ele_isPF',ele_type='Both',isPF=ak.values_astype(fake_both.isPF,int),weight=1)
        h.fill('fake_ele_numTrkHits',ele_type='Both',numHits=fake_both.numTrackerHits,weight=1)
        h.fill('fake_ele_numPixHits',ele_type='Both',numHits=fake_both.numPixHits,weight=1)
        h.fill('fake_ele_numStripHits',ele_type='Both',numHits=fake_both.numStripHits,weight=1)
        h.fill("fake_ele_sigmaIetaIeta",ele_type='Both',sieie=fake_both.full55sigmaIetaIeta,weight=1)
        h.fill("fake_ele_absdEtaSeed",ele_type='Both',detaSeed=fake_both.absdEtaSeed,weight=1)
        h.fill("fake_ele_absdPhiIn",ele_type='Both',abs_dphi=np.abs(fake_both.absdPhiIn),weight=1)
        h.fill("fake_ele_HOverE",ele_type='Both',hoe=fake_both.HoverE,weight=1)
        h.fill("fake_ele_1E1p",ele_type='Both',emp=fake_both.abs1overEm1overP,weight=1)
        h.fill("fake_ele_expMissing",ele_type='Both',missing=fake_both.expMissingInnerHits,weight=1)
        h.fill("fake_ele_passConvVeto",ele_type='Both',passVeto=fake_both.conversionVeto,weight=1)
        h.fill("fake_ele_IDscore",ele_type='Both',id=fake_both.IDscore,weight=1)

        h.fill("fake_ele_pt",ele_type='L',pt=fake_lpt.pt,weight=1)
        h.fill("fake_ele_eta",ele_type='L',eta=fake_lpt.eta,weight=1)
        h.fill("fake_ele_phi",ele_type='L',phi=fake_lpt.phi,weight=1)
        h.fill("fake_ele_dxy",ele_type='L',dxy=fake_lpt.dxy,weight=1)
        h.fill("fake_ele_dxySignif",ele_type='L',signif=np.abs(fake_lpt.dxy)/fake_lpt.dxyErr,weight=1)
        h.fill("fake_ele_trkChi2",ele_type='L',chi2=fake_lpt.trkChi2,weight=1)
        h.fill("fake_ele_trkProb",ele_type='L',prob=fake_lpt.trkProb)
        h.fill("fake_ele_trkRelIso",ele_type='L',iso=fake_lpt.trkRelIso,weight=1)
        h.fill("fake_ele_calRelIso",ele_type='L',iso=fake_lpt.calRelIso,weight=1)
        h.fill("fake_ele_PFRelIso",ele_type='L',iso=fake_lpt.PFRelIso,weight=1)
        h.fill("fake_ele_miniRelIso",ele_type='L',iso=fake_lpt.miniRelIso,weight=1)
        h.fill("fake_ele_mindRJets",ele_type='L',dR=fake_lpt.mindRj,weight=1)
        h.fill("fake_ele_mindPhiJets",ele_type='L',abs_dphi=np.abs(fake_lpt.mindPhiJ),weight=1)
        h.fill('fake_ele_isPF',ele_type='L',isPF=ak.values_astype(fake_lpt.isPF,int),weight=1)
        h.fill('fake_ele_numTrkHits',ele_type='L',numHits=fake_lpt.numTrackerHits,weight=1)
        h.fill('fake_ele_numPixHits',ele_type='L',numHits=fake_lpt.numPixHits,weight=1)
        h.fill('fake_ele_numStripHits',ele_type='L',numHits=fake_lpt.numStripHits,weight=1)
        h.fill("fake_ele_sigmaIetaIeta",ele_type='L',sieie=fake_lpt.full55sigmaIetaIeta,weight=1)
        h.fill("fake_ele_absdEtaSeed",ele_type='L',detaSeed=fake_lpt.absdEtaSeed,weight=1)
        h.fill("fake_ele_absdPhiIn",ele_type='L',abs_dphi=np.abs(fake_lpt.absdPhiIn),weight=1)
        h.fill("fake_ele_HOverE",ele_type='L',hoe=fake_lpt.HoverE,weight=1)
        h.fill("fake_ele_1E1p",ele_type='L',emp=fake_lpt.abs1overEm1overP,weight=1)
        h.fill("fake_ele_expMissing",ele_type='L',missing=fake_lpt.expMissingInnerHits,weight=1)
        h.fill("fake_ele_passConvVeto",ele_type='L',passVeto=fake_lpt.conversionVeto,weight=1)
        h.fill("fake_ele_IDscore",ele_type='L',id=fake_lpt.IDscore,weight=1)

        #
        h.fill("match_vtx_dR",vtype=match_vtx.typ,dR=match_vtx.dR,weight=1)
        h.fill("match_vtx_mindxy",vtype=match_vtx.typ,dxy=match_vtx.min_dxy,weight=1)
        h.fill("match_vtx_vxy1",vtype=match_vtx.typ,vxy=match_vtx.vxy,weight=1)
        h.fill("match_vtx_vxy10",vtype=match_vtx.typ,vxy=match_vtx.vxy,weight=1)
        h.fill("match_vtx_vxy100",vtype=match_vtx.typ,vxy=match_vtx.vxy,weight=1)
        h.fill("match_vtx_leadpT",vtype=match_vtx.typ,pt=np.maximum(match_vtx.e1.pt,match_vtx.e2.pt),weight=1)
        h.fill("match_vtx_METdPhi",vtype=match_vtx.typ,abs_dphi=np.abs(match_vtx.METdPhi),weight=1)
        h.fill('match_vtx_mindRj',vtype=match_vtx.typ,dR=match_vtx.mindRj,weight=1)
        h.fill('match_vtx_chi2',vtype=match_vtx.typ,chi2=match_vtx.reduced_chi2,weight=1)
        h.fill('match_vtx_mass',vtype=match_vtx.typ,mass=match_vtx.m,weight=1)
        h.fill('match_vtx_mindPhiJ',vtype=match_vtx.typ,abs_dphi=np.abs(match_vtx.mindPhiJ),weight=1)
        h.fill('match_vtx_sign',vtype=match_vtx.typ,sign=match_vtx.sign,weight=1)
        h.fill('match_vtx_pt',vtype=match_vtx.typ,pt=match_vtx.pt,weight=1)
        h.fill('match_vtx_eta',vtype=match_vtx.typ,eta=match_vtx.eta,weight=1)
        h.fill('match_vtx_phi',vtype=match_vtx.typ,phi=match_vtx.phi,weight=1)
        h.fill('match_vtx_type',vtype=match_vtx.typ,weight=1)
        h.fill("match_vtx_minEleDrJ",vtype=match_vtx.typ,dR=np.minimum(match_vtx.e1.mindRj,match_vtx.e2.mindRj),weight=1)
        h.fill("match_vtx_minEleDPhiJ",vtype=match_vtx.typ,abs_dphi=np.minimum(match_vtx.e1.mindPhiJ,match_vtx.e2.mindPhiJ),weight=1)
        h.fill("match_vtx_mass_low",vtype=match_vtx.typ,mass_low=match_vtx.m,weight=1)
        h.fill("match_vtx_mindxy_low",vtype=match_vtx.typ,mindxy_low=match_vtx.min_dxy,weight=1)
        h.fill("match_vtx_sign_etaProd",vtype=match_vtx.typ,sign_etaProd=ak.values_astype(np.sign(match_vtx.e1.eta*match_vtx.e2.eta),int),weight=1)
        h.fill("match_vtx_CosThetaColl",vtype=match_vtx.typ,cosTheta=match_vtx_cosTheta,weight=1)
        h.fill("match_vtx_LxyCosThetaColl",vtype=match_vtx.typ,LxyCosTheta=match_vtx_cosTheta*match_vtx.vxy,weight=1)
        h.fill("match_vtx_LxyCosThetaCollZoom",vtype=match_vtx.typ,LxyCosThetaZoom=match_vtx_cosTheta*match_vtx.vxy,weight=1)
        h.fill("match_vtx_LxyCosThetaCollZoomZoom",vtype=match_vtx.typ,LxyCosThetaZoomZoom=match_vtx_cosTheta*match_vtx.vxy,weight=1)
        h.fill("match_vtx_eleDphi",vtype=match_vtx.typ,abs_dphi=match_vtx.eleDphi,weight=1)
        h.fill("match_vtx_chi2Rank",vtype=match_vtx.typ,chi2Rank=match_vtx_chi2Rank,weight=1)

        #
        h.fill("match_vtx_gen_dR",vtype=match_vtx.typ,dR=match_vtx_genObj.dr,weight=1)
        h.fill("match_vtx_gen_vxy1",vtype=match_vtx.typ,vxy=match_vtx_genObj.vxy,weight=1)
        h.fill("match_vtx_gen_vxy10",vtype=match_vtx.typ,vxy=match_vtx_genObj.vxy,weight=1)
        h.fill("match_vtx_gen_vxy100",vtype=match_vtx.typ,vxy=match_vtx_genObj.vxy,weight=1)
        h.fill("match_vtx_gen_METdPhi",vtype=match_vtx.typ,abs_dphi=np.abs(match_vtx_genObj.METdPhi),weight=1)
        h.fill("match_vtx_gen_mass",vtype=match_vtx.typ,mass=match_vtx_genObj.mass,weight=1)
        h.fill("match_vtx_gen_pt",vtype=match_vtx.typ,pt=match_vtx_genObj.pt,weight=1)
        h.fill("match_vtx_gen_eta",vtype=match_vtx.typ,eta=match_vtx_genObj.eta,weight=1)
        h.fill("match_vtx_gen_phi",vtype=match_vtx.typ,phi=match_vtx_genObj.phi,weight=1)
        
        #
        h.fill("fake_vtx_dR",vtype=fake_vtx.typ,dR=fake_vtx.dR,weight=1)
        h.fill("fake_vtx_mindxy",vtype=fake_vtx.typ,dxy=fake_vtx.min_dxy,weight=1)
        h.fill("fake_vtx_vxy1",vtype=fake_vtx.typ,vxy=fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_vxy10",vtype=fake_vtx.typ,vxy=fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_vxy100",vtype=fake_vtx.typ,vxy=fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_leadpT",vtype=fake_vtx.typ,pt=np.maximum(fake_vtx.e1.pt,fake_vtx.e2.pt),weight=1)
        h.fill("fake_vtx_METdPhi",vtype=fake_vtx.typ,abs_dphi=np.abs(fake_vtx.METdPhi),weight=1)
        h.fill('fake_vtx_mindRj',vtype=fake_vtx.typ,dR=fake_vtx.mindRj,weight=1)
        h.fill('fake_vtx_chi2',vtype=fake_vtx.typ,chi2=fake_vtx.reduced_chi2,weight=1)
        h.fill('fake_vtx_mass',vtype=fake_vtx.typ,mass=fake_vtx.m,weight=1)
        h.fill('fake_vtx_mindPhiJ',vtype=fake_vtx.typ,abs_dphi=np.abs(fake_vtx.mindPhiJ),weight=1)
        h.fill('fake_vtx_sign',vtype=fake_vtx.typ,sign=fake_vtx.sign,weight=1)
        h.fill('fake_vtx_pt',vtype=fake_vtx.typ,pt=fake_vtx.pt,weight=1)
        h.fill('fake_vtx_eta',vtype=fake_vtx.typ,eta=fake_vtx.eta,weight=1)
        h.fill('fake_vtx_phi',vtype=fake_vtx.typ,phi=fake_vtx.phi,weight=1)
        h.fill('fake_vtx_type',vtype=fake_vtx.typ,weight=1)
        h.fill("fake_vtx_minEleDrJ",vtype=fake_vtx.typ,dR=np.minimum(fake_vtx.e1.mindRj,fake_vtx.e2.mindRj),weight=1)
        h.fill("fake_vtx_minEleDPhiJ",vtype=fake_vtx.typ,abs_dphi=np.minimum(fake_vtx.e1.mindPhiJ,fake_vtx.e2.mindPhiJ),weight=1)
        h.fill("fake_vtx_mass_low",vtype=fake_vtx.typ,mass_low=fake_vtx.m,weight=1)
        h.fill("fake_vtx_mindxy_low",vtype=fake_vtx.typ,mindxy_low=fake_vtx.min_dxy,weight=1)
        h.fill("fake_vtx_sign_etaProd",vtype=fake_vtx.typ,sign_etaProd=ak.values_astype(np.sign(fake_vtx.e1.eta*fake_vtx.e2.eta),int),weight=1)
        h.fill("fake_vtx_CosThetaColl",vtype=fake_vtx.typ,cosTheta=fake_vtx_cosTheta,weight=1)
        h.fill("fake_vtx_LxyCosThetaColl",vtype=fake_vtx.typ,LxyCosTheta=fake_vtx_cosTheta*fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_LxyCosThetaCollZoom",vtype=fake_vtx.typ,LxyCosThetaZoom=fake_vtx_cosTheta*fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_LxyCosThetaCollZoomZoom",vtype=fake_vtx.typ,LxyCosThetaZoomZoom=fake_vtx_cosTheta*fake_vtx.vxy,weight=1)
        h.fill("fake_vtx_eleDphi",vtype=fake_vtx.typ,abs_dphi=fake_vtx.eleDphi,weight=1)

        #
        h.fill("PFMET",met=events.PFMET.pt,weight=wgt)
        h.fill("PFMET_vs_genMET",met=events.PFMET.pt,genmet=events.GenMET.pt,weight=1)
        h.fill("jetMETdPhi",abs_dphi=np.abs(events.PFJet.METdPhi[:,0]),weight=wgt)
        h.fill("minJetMETdPhi",abs_dphi=ak.min(np.abs(events.PFJet.METdPhi),axis=1),weight=wgt)
        h.fill("nJets",nJets=ak.count(events.PFJet.pt,axis=1),weight=wgt)
        h.fill("jetMETratio",jetMETratio=events.PFJet.pt[:,0]/events.PFMET.pt,weight=wgt)

        sel_vtx = events.sel_vtx
        h.fill('sel_vtx_isMatched',matched=sel_vtx.isMatched,weight=wgt)
        h.fill('sel_vtx_purity',purity='total',weight=np.sum(wgt))
        h.fill('sel_vtx_purity',purity='n_matched',weight=np.sum(wgt[sel_vtx.isMatched]))
        h.fill('sel_vtx_purity',purity='n_eeReco',weight=np.sum(wgt[events.signalReconstructed]))
        h.fill('sel_vtx_purity',purity='n_vtxReco',weight=np.sum(wgt[events.signalReconstructed & ak.any(events.good_vtx.isMatched,axis=1)]))
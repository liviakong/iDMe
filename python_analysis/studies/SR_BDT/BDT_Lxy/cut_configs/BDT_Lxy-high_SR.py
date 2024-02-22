import numpy as np
import awkward as ak
import sys
sys.path.append("../../analysisTools/")
import analysisSubroutines as routines

def cut5(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut5"
    desc = "No b-tagged jets"
    plots = False
    bTag = events.PFJet.bTag
    # DeepFlavour working points for UL samples
    if info["year"] == 2018:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        wp = 0.0490 # loose
        #wp = 0.2783 # medium
        #wp = 0.7100 # tight
    if info["year"] == 2017:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
        wp = 0.0532 # loose
        #wp = 0.3040 # medium
        #wp = 0.7476 # tight
    if info["year"] == "2016_preVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
        wp = 0.0508 # loose
        #wp = 0.2598 # medium
        #wp = 0.6502 # tight
    if info["year"] == "2016_postVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
        wp = 0.0480 # loose
        #wp = 0.2489 # medium
        #wp = 0.6377 # tight
    pass_bTag = bTag > wp
    n_bTag_Jets = ak.count(events.PFJet[pass_bTag].pt,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = False
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "LxyBin"
    desc = "Lxy > 0.5"
    plots = True
    cut = events.sel_vtx.vxy > 0.5
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "BDT"
    plots = True

    variables = ['sel_vtx_sign', 'sel_vtx_chi2','sel_vtx_METdPhi','sel_vtx_m','sel_vtx_dR','sel_vtx_minDxy','vxy_signif'] # BDT variables

    vxy_thres = 0.5
    vxy_mask = events["sel_vtx","vxy"] > vxy_thres

    thres_high_vxy = 0.9524
    thres_low_vxy = 0.8267
    
    thres_arr = ak.where(vxy_mask, thres_high_vxy, thres_low_vxy)

    if len(events) != 0:
        input = routines.makeBDTv2Inputs(events)

        model_high_vxy = './models/BDT_highLxy_thres_0p5.json'
        model_low_vxy = './models/BDT_lowLxy_thres_0p5.json'
        
        score_BDT_high_vxy = routines.getBDTscore(input, model_high_vxy)
        score_BDT_low_vxy = routines.getBDTscore(input, model_low_vxy)

        #score_BDT = ak.fill_none(ak.mask(score_BDT_high_vxy, vxy_mask), 0) + ak.fill_none(ak.mask(score_BDT_low_vxy, ~vxy_mask), 0)
        #cut = score_BDT > thres_arr

        cut = vxy_mask & (score_BDT_high_vxy > thres_high_vxy)

        print(f'Pass: {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []
    
    return events[cut], name, desc, plots
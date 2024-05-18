import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Preselection"
    plots = True
    return events, name, desc, plots

def cut1(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut1"
    desc = "No b-tagged jets"
    plots = True
    bTag = events.PFJet.bTag
    # DeepFlavour working points for UL samples
    if info["year"] == 2018:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        #wp = 0.0490 # loose
        wp = 0.2783 # medium
        #wp = 0.7100 # tight
    if info["year"] == 2017:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
        #wp = 0.0532 # loose
        wp = 0.3040 # medium
        #wp = 0.7476 # tight
    if info["year"] == "2016_preVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
        #wp = 0.0508 # loose
        wp = 0.2598 # medium
        #wp = 0.6502 # tight
    if info["year"] == "2016_postVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
        #wp = 0.0480 # loose
        wp = 0.2489 # medium
        #wp = 0.6377 # tight
    pass_bTag = bTag > wp
    n_bTag_Jets = ak.count(events.PFJet[pass_bTag].pt,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "Leading jet |eta| < 2.4"
    plots = True
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

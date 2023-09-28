import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Has valid vertex"
    plots = False
    return events, name, desc, plots

def cut1(events,info):
    name = "cut1"
    desc = "Pass MET Filters"
    plots = False
    cut = events.METFiltersFailBits == 0
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "HEM Veto"
    plots = False
    if info["year"] == 2018:
        cut = ~events.HEM.flag
        return events[cut], name, desc, plots
    else:
        return events, name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "MET Trigger (120 GeV)"
    plots = False
    if info["year"] == 2016:
        cut = (events.trigFired16 & (1<<9)) == (1<<9)
    if info["year"] == 2017:
        cut = (events.trigFired17 & (1<<9)) == (1<<9)
    if info["year"] == 2018:
        cut = (events.trigFired18 & (1<<13)) == (1<<13)
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "MET > 200 GeV"
    plots = False
    cut = events.PFMET.pt > 200
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "0 < nJets < 3 (pT > 30 GeV)"
    plots = False
    nJets = ak.count(events.PFJet.pt,axis=1)
    cut = (nJets > 0) & (nJets < 3)
    return events[cut], name, desc, plots

def cut6(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut6"
    desc = "No b-tagged jets"
    plots = False
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

def cut7(events,info):
    name = "cut7"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "Leading jet pT > 80 GeV"
    plots = False
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = False
    cut = events.PFJet.METdPhi[:,0] > 1.5
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = False
    cut = ak.all(events.PFJet.METdPhi > 0.75,axis=1)
    return events[cut], name, desc, plots
    
def cut11(events,info):
    name = "cut11"
    desc = "OSSF"
    plots = True
    cut = events.sel_vtx.sign == -1
    return events[cut], name, desc, plots

def cut12(events,info):
    name = "cut12"
    desc = "Vertex chi2/df < 5"
    plots = True
    cut = events.sel_vtx.reduced_chi2 < 5
    return events[cut], name, desc, plots

def cut13(events,info):
    name = "cut13"
    desc = "dPhi(MET,vtx) < 2.5"
    plots = True
    cut = events.sel_vtx.METdPhi < 2.5
    return events[cut], name, desc, plots

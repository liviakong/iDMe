import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Preselection"
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
    cut = (info['year'] == 2018) & (~events.HEM.flag)
    return events[cut], name, desc, plots

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
    plots = True
    cut = events.PFMET.correctedPt > 200
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "|Calo MET - PF MET|/Calo MET < 1.0"
    plots = False
    cut = (np.abs(events.CaloMET.pt - events.PFMET.correctedPt) / events.CaloMET.pt) < 1.0
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "0 < nJets < 3 (pT > 30 GeV)"
    plots = False
    nJets = ak.count(events.PFJet.corrPt,axis=1)
    cut = (nJets > 0) & (nJets < 3)
    return events[cut], name, desc, plots

def cut7(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut7"
    desc = "No b-tagged jets"
    plots = False
    bTag = events.PFJet.bTag
    if info["year"] == 2018:
        wp_med = 0.4168
    if info["year"] == 2017:
        wp_med = 0.4506
    if info["year"] == "2016_preVFP":
        wp_med = 0.6001
    if info["year"] == "2016_postVFP":
        wp_med = 0.5847
    pass_bTag = bTag > wp_med
    n_bTag_Jets = ak.count(events.PFJet[pass_bTag].pt,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.corrEta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.corrPt[:,0] > 80
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = False
    cut = events.JetMETdPhi[:,0] > 1.5
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(events.JetMETdPhi > 0.75,axis=1)
    return events[cut], name, desc, plots
    
def cut12(events,info):
    name = "cut12"
    desc = "OSSF"
    plots = True
    cut = (events.sel_e1.charge * events.sel_e2.charge) == -1
    return events[cut], name, desc, plots

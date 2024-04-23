import numpy as np
import awkward as ak

def cut1(events,info):
    name = "cut1"
    desc = "1+ good vertices"
    plots = True
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut1"
    desc = "0 < nJets < 3"
    plots = True
    nJets = ak.count(events.PFJet.pt,axis=1)
    cut = (nJets > 0) & (nJets < 3)
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "MET Filters"
    plots = True
    cut = events.METFiltersFailBits == 0
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "HEM Veto"
    plots = True
    if info['year'] == 2018:
        cut = ~events.HEM.flag
        return events[cut], name, desc, plots
    else:
        return events, name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "MET Trigger 120"
    plots = True
    cut = events.trig.HLT_PFMET120_PFMHT120_IDTight == 1
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "MET > 200"
    plots = True
    cut = events.PFMET.pt > 200
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "No b-tagged (med)"
    plots = True
    n_bTag_Jets = ak.count_nonzero(events.PFJet.passMedID,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "Leading jet |eta| < 2.4"
    plots = True
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

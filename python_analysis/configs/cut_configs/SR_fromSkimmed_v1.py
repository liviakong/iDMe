import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Leading jet |eta| < 2.4"
    plots = True
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut1(events,info):
    name = "cut1"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots
    
def cut4(events,info):
    name = "cut4"
    desc = "OSSF"
    plots = True
    cut = events.sel_vtx.sign == -1
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "Vertex chi2/df < 5"
    plots = True
    cut = events.sel_vtx.reduced_chi2 < 5
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "dPhi(MET,vtx) < 2.5"
    plots = True
    cut = np.abs(events.sel_vtx.METdPhi) < 2.5
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "Vertex chi2/df < 3"
    plots = True
    cut = events.sel_vtx.reduced_chi2 < 3
    return events[cut], name, desc, plots
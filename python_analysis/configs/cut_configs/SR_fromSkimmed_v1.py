import numpy as np
import awkward as ak

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
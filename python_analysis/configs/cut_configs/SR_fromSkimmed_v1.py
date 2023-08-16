import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut1(events,info):
    name = "cut1"
    desc = "Leading jet pT > 80 GeV"
    plots = False
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = False
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = False
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots
    
def cut4(events,info):
    name = "cut4"
    desc = "OSSF"
    plots = False
    cut = events.sel_vtx.sign == -1
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "Vertex chi2/df < 5"
    plots = False
    cut = events.sel_vtx.reduced_chi2 < 5
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "dPhi(MET,vtx) < 2.5"
    plots = False
    cut = np.abs(events.sel_vtx.METdPhi) < 2.5
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "Vertex chi2/df < 3"
    plots = True
    cut = events.sel_vtx.reduced_chi2 < 3
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "dPhi(MET,vtx) < 1.0"
    plots = True
    cut = np.abs(events.sel_vtx.METdPhi) < 1.0
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "m(ee) < 20 GeV"
    plots = True
    cut = events.sel_vtx.m < 20
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "dR(ee) < 1.0"
    plots = True
    cut = events.sel_vtx.dR < 1.0
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "min(dxy) > 0.01"
    plots = True
    min_dxy = np.minimum(np.abs(events.sel_vtx.e1.dxy),np.abs(events.sel_vtx.e2.dxy))
    cut = min_dxy > 0.01
    return events[cut], name, desc, plots
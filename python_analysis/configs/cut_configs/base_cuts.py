import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    info = "No cuts"
    return events, name, info

def cut1(events,info):
    name = "cut1"
    info = "Pass MET Filters"
    cut = events.METFiltersFailBits == 0
    return events[cut], name, info

def cut2(events,info):
    name = "cut2"
    info = "HEM Veto"
    cut = (info['year'] == 2018) & (~events.HEM.flag)
    return events[cut], name, info

def cut3(events,info):
    name = "cut3"
    info = "MET Trigger (120 GeV)"
    if info["year"] == 2016:
        cut = (events.trigFired16 & (1<<9)) == (1<<9)
    if info["year"] == 2017:
        cut = (events.trigFired17 & (1<<9)) == (1<<9)
    if info["year"] == 2018:
        cut = (events.trigFired18 & (1<<13)) == (1<<13)
    return events[cut], name, info

def cut4(events,info):
    name = "cut4"
    info = "MET > 200 GeV"
    cut = events.PFMET.correctedPt > 200
    return events[cut], name, info

def cut5(events,info):
    name = "cut5"
    info = "|Calo MET - PF MET|/Calo MET < 1.0"
    cut = (np.abs(events.CaloMET.pt - events.PFMET.correctedPt) / events.CaloMET.pt) < 1.0
    return events[cut], name, info
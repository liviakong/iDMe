import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Dummy cut"
    plots = True
    return events, name, desc, plots

def cut1(events,info):
    name = 'cut1'
    desc = 'goodVtx default'
    plots = True
    info['defineGoodVertices'](events,version='default',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut2(events,info):
    name = 'cut2'
    desc = 'goodVtx v1'
    plots = True
    info['defineGoodVertices'](events,version='v1',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut3(events,info):
    name = 'cut3'
    desc = 'goodVtx v2'
    plots = True
    info['defineGoodVertices'](events,version='v2',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut4(events,info):
    name = 'cut4'
    desc = 'goodVtx v3'
    plots = True
    info['defineGoodVertices'](events,version='v3',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut5(events,info):
    name = 'cut5'
    desc = 'goodVtx v4'
    plots = True
    info['defineGoodVertices'](events,version='v4',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots

def cut6(events,info):
    name = 'cut6'
    desc = 'goodVtx v5'
    plots = True
    info['defineGoodVertices'](events,version='v5',ele_id='basic')
    info['selectBestVertex'](events)
    cut = events.nGoodVtx > 0
    return events[cut], name, desc, plots
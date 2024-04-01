from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        # quantities associated w/ selected vertex
        "bdtScore" : Hist(samp,cut,bdtScore,storage=hist.storage.Weight()),
        "gen_met_noWgt" :  Hist(samp,cut,met,storage=hist.storage.Weight()),
        "gen_met" :  Hist(samp,cut,met,storage=hist.storage.Weight()),
        "gen_leadjet_pt" : Hist(samp,cut,jet_pt,storage=hist.storage.Weight()),
        "gen_leadjet_pt_noWgt" : Hist(samp,cut,jet_pt,storage=hist.storage.Weight())

    }
    return histograms

subroutines = []

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    wgt = events.eventWgt/sum_wgt
    if 'BDTScore' in events.fields:
        histos["bdtScore"].fill(samp=samp,cut=cut,score=events.BDTScore,weight=wgt)
    if info['type'] == "signal":
        histos['gen_met_noWgt'].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=1)
        histos['gen_met'].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=wgt)
        histos['gen_leadjet_pt'].fill(samp=samp,cut=cut,pt=events.GenJet.pt[:,0],weight=wgt)
        histos['gen_leadjet_pt_noWgt'].fill(samp=samp,cut=cut,pt=events.GenJet.pt[:,0],weight=1)



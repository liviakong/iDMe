from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        # quantities associated w/ selected vertex
        "bdtScore" : Hist(samp,cut,bdtScore,storage=hist.storage.Weight()),
    }
    return histograms

subroutines = []

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    wgt = events.eventWgt/sum_wgt
    if 'BDTScore' in events.fields:
        histos["bdtScore"].fill(samp=samp,cut=cut,score=events.BDTScore,weight=wgt)

import numpy as np
import awkward as ak
from myHisto import myHisto

def make_histograms(info):
    h = myHisto()
        
    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    return
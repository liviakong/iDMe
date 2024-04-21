import hist
import coffea.util as util
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re
import os
import sys
import mplhep as hep

def makeNMinus1(h1,h2,lessThan=False):
    assert len(h1.axes)==1 and h1.axes == h2.axes
    ax = h1.axes[0]
    if lessThan:
        s1 = np.cumsum(h1.counts(flow=True))
        s2 = np.cumsum(h2.counts(flow=True))
    else:
        s1 = np.cumsum(h1.counts(flow=True)[::-1])[::-1]
        s2 = np.cumsum(h2.counts(flow=True)[::-1])[::-1]
    x = ax.edges
    dx = x[1]-x[0]
    x = np.append(x,[x[-1]+dx])
    signif = np.where(s2>0,s1/np.sqrt(s2),-1)
    signif[(signif==-1) & (s1>0)] = np.inf
    signif[(signif==-1) & (s1==0)] = 0
    return x, signif

def makeNMinus1_multiBkg(hnum,hdens,lessThan=False):
    assert len(hnum.axes)==1 and hnum.axes == hdens[0].axes
    ax = hnum.axes[0]
    if lessThan:
        snum = np.cumsum(hnum.counts(flow=True))
        sden = sum([np.cumsum(hd.counts(flow=True)) for hd in hdens])
    else:
        snum = np.cumsum(hnum.counts(flow=True)[::-1])[::-1]
        sden = sum([np.cumsum(hd.counts(flow=True)[::-1])[::-1] for hd in hdens])
    x = ax.edges
    dx = x[1]-x[0]
    x = np.append(x,[x[-1]+dx])
    signif = np.where(sden>0,snum/np.sqrt(sden),-1)
    signif[(signif==-1) & (snum>0)] = np.inf
    signif[(signif==-1) & (snum==0)] = 0
    return x, signif

def makeCutEff(h,lessThan=False):
    ax = h.axes[0]
    if lessThan:
        s = np.cumsum(h.counts(flow=True)[:-1])
    else:
        s = np.cumsum(h.counts(flow=True)[1:][::-1])[::-1]
    x = ax.edges
    #dx = x[1]-x[0]
    #x = np.append(x,[x[-1]+dx])
    eff = s/np.sum(h.counts(flow=True))
    return x, eff, s
    

def getSampleInfo(histos,hname="ele_kinematics"):
    samps = [s.name for s in histos[hname].axis("sample").identifiers()]
    info = {}
    masses = []
    cts = []
    for s in samps:
        ct = re.findall("ctau-(\d+)",s)[0]
        m, dm = re.findall("Mchi-(\d+p\d)_dMchi-(\d+p\d)",s)[0]
        m = m.replace("p",".")
        dm = dm.replace("p",".")
        entry = "{0}-{1}".format(m,dm)
        
        if entry not in info.keys():
            info[entry] = {}
        info[entry][ct] = s

        if entry not in masses:
            masses.append(entry)
        if ct not in cts:
            cts.append(ct)
    return info, masses, cts

def reduceSampleName(name,lifetime=False,mass=False,full=False,verbosity=0):
    output = name
    m, dm = re.findall("Mchi-(\d+p\d)_dMchi-(\d+p\d)",name)[0]
    ct = re.findall("ctau-(\d+)",name)[0]
    m = m.replace("p",".")
    dm = dm.replace("p",".")
    if full:
        output = r"$m_\chi = {0}$ GeV, $\Delta m_\chi = {1}$ GeV, $c\tau = {2}$ mm".format(m,dm,ct)
    if lifetime:
        output = r"$c\tau = {0}$ mm".format(ct)
    if mass:
        if verbosity == 0: output = r"$({0}, {1})$ GeV".format(m,dm)
        if verbosity == 1: output = r"$(m_\chi, \Delta m_\chi) = ({0}, {1})$ GeV".format(m,dm)
        if verbosity == 2: output = r"$m_\chi = {0}$ GeV, $\Delta m_\chi = {1}$ GeV".format(m,dm)
    return output

def signalPoint(name):
    a = re.search('Mchi-(\d+p\d+)_dMchi-(\d+p\d+)_ctau-(\d+)',name)
    mchi = float(a.group(1).replace("p","."))
    dmchi = float(a.group(2).replace("p","."))
    ctau = float(a.group(3).replace("p","."))
    m1 = mchi - dmchi/2
    m2 = mchi + dmchi/2
    delta = dmchi/m1
    return {"mchi":mchi, "dmchi":dmchi, "ctau":ctau, "m1":m1, "m2":m2, "delta":delta, "name":name}

def getCut(label,n=2):
    name = ""
    while label[0:n]!=label[n:2*n] and n<len(label):
        name=label[0:n+1]
        n+=1
    return name

def getHTlow(sampName):
    ht = int(re.search("HT(\d+)to",sampName).group(1))
    return ht

def setDefaultStyle(fontsize=14):
    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["figure.figsize"] = (10,8)

def hget(h,samp,cut):
    return h[{"samp":samp,"cut":cut}]

def makeCDF(h,start,stop,bins=100,right=True,nevents=False):
    x = np.linspace(start,stop,bins)
    n_tot = h.sum(flow=True).value
    if right:
        yields = np.array([h[complex(f"{xi}j")::sum].value for xi in x])
    else:
        yields = np.array([h[:complex(f"{xi}j"):sum].value for xi in x])
    effs = yields/n_tot
    if nevents:
        return x,yields
    else:
        return x, effs

def overlay(h,overlay,label_key=None,**kwargs):
    axes = h.axes
    targ = None
    for a in axes:
        if a.name == overlay:
            targ = a
    if targ is None:
        print("can't find overlay axis!")
        return
    n_overlay = len(targ.centers)
    labels = [targ.value(i) for i in range(n_overlay)]
    histos = [h[{overlay:l}] for l in labels]
    if label_key is not None:
        labels = [label_key[targ.value(i)] for i in range(n_overlay)]
    hep.histplot(histos,label=labels,**kwargs)

from mplhep.styles.cms import cmap_petroff
bkg_cmap = {
    "QCD":cmap_petroff[0],
    "WJets":cmap_petroff[1],
    "ZJets":cmap_petroff[2],
    "DY":cmap_petroff[3],
    "Top":cmap_petroff[4],
    "Multiboson":cmap_petroff[5]
}

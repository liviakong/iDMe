import coffea.hist as hist
import coffea.util as util
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re
import os
import sys

def loadHistoFiles(location):
    files = [f for f in os.listdir(location) if ".coffea" in f]
    histos = {}
    for f in files:
        htemp = util.load(location+"/"+f)[0]
        for hname in list(htemp.keys()):
            if hname not in histos.keys():
                histos[hname] = htemp[hname].copy()
            else:
                histos[hname] += htemp[hname]
    return histos

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

def setDefaultStyle():
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["figure.figsize"] = (10,8)
from __future__ import with_statement
import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from mySchema import MySchema
from coffea import processor
import coffea.hist as hist
from coffea.nanoevents.methods import vector
import uproot
import awkward as ak
ak.behavior.update(vector.behavior)
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import time
import importlib
import pandas as pd
from XRootD import client
import re
NanoAODSchema.warn_missing_crossrefs = False
import analysisSubroutines as routines
import sys

match_names = {"Default":"match0","lowpt":"match1"}
vxy_range = {1:[0,20],10:[0,50],100:[0,50],1000:[0,50]}
vxy_rebin = {1:5,10:20,100:20,1000:20}

class Analyzer:
    def __init__(self,fileList,histoList,cuts,max_samples=-1,max_files_per_samp=-1):
        # load in file config
        if type(fileList) == str and ".json" in fileList:
            with open(fileList) as f:
                self.fileList = json.load(f)
        else:
            self.fileList = fileList
        
        #load in histogram config
        if "/" in histoList: # if cut file is in a different directory
            sys.path.append("/".join(histoList.split("/")[:-1]))
            self.histoFile = histoList.split("/")[-1].split(".")[0]
        else: # cut file is in the same directory (e.g. running on condor)
            self.histoFile = histoList.split(".")[0]

        # load in cuts config (should be a path to a .py file with cut methods)
        self.cuts = cuts

        self.sample_names = [] # list of sample names, readable names are generated from data in the fileList json
        self.sample_locs = {} # dictionary mapping sample name to file/directory location
        self.sample_info = {} # dictionary with sample metadata
        self.max_samples = max_samples
        self.max_files_per_samp = max_files_per_samp
        self.nEvents = {}
        self.nEventsProcessed = {}
        self.totalEvents = 0
        self.mode = None
        
        self.loadFiles()
    
    def loadFiles(self):
        loaded = 0
        for sample in self.fileList:
            if self.max_samples > 0 and loaded == self.max_samples:
                break
            mode = sample['type']
            if mode == 'signal':
                name = "sig_{0}_Mchi-{1}_dMchi-{2}_ctau-{3}".format(sample['year'],str(sample['Mchi']).replace(".","p"),str(sample['dMchi']).replace(".","p"),sample['ctau'])
            elif mode == 'bkg':
                name = "bkg_{0}_{1}".format(sample['year'],sample['name'])
            elif mode == 'data':
                name = "data_{0}_{1}".format(sample['year'],sample['name'])
            
            if self.mode is None:
                self.mode = mode
            else:
                if self.mode != mode:
                    print("Error! You're mixing samples of differing types (e.g. signal and bkg, signal and data, etc)")
                    print("Please split up different kinds of samples into different configs")
                    exit()
            
            loc = sample['location']
            if '.root' in loc:
                # if the location is just a single file, load it in
                self.sample_locs[name] = [sample['location']]
            elif 'fileset' in sample.keys():
                self.sample_locs[name] = sample['fileset']
            else:
                # if the location is a directory, use the xrootd client to get a list of files
                xrdClient = client.FileSystem("root://cmseos.fnal.gov")
                status, flist = xrdClient.dirlist(loc)
                fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if '.root' in item.name]
                if self.max_files_per_samp > 0:
                    fullList = fullList[:self.max_files_per_samp] if len(fullList) > self.max_files_per_samp else fullList
                self.sample_locs[name] = fullList
            
            self.sample_info[name] = sample
            self.sample_names.append(name)
            loaded += 1

    def process(self,treename='ntuples/outT',execr="iterative"):
        fileset = self.sample_locs
        proc = iDMeProcessor(self.sample_names,self.sample_info,self.sample_locs,self.histoFile,self.cuts,mode=self.mode)
        if execr == "iterative":
            executor = processor.iterative_executor
            executor_args = {"schema":MySchema}
        elif execr == "futures":
            executor = processor.futures_executor
            executor_args = {"workers": 4,
                        "savemetrics": True,
                        "schema": MySchema,
                        "align_clusters": True
                        }
        accumulator = processor.run_uproot_job(fileset,
                              treename=treename,
                              processor_instance=proc,
                              executor=executor,
                              executor_args=executor_args)
        return accumulator

class iDMeProcessor(processor.ProcessorABC):
    def __init__(self,samples,sampleInfo,fileSet,histoFile,cutFile,mode='signal'):
        self.samples = samples
        self.sampleInfo = sampleInfo
        self.sampleLocs = fileSet
        self.mode = mode
        
        # load in histogram config
        histoMod = importlib.import_module(histoFile)
        histos = histoMod.histograms
        self.histoFill = histoMod.fillHistos
        self.subroutines = histoMod.subroutines

        # make a cutflow dictionary for the output
        cutflows = processor.dict_accumulator({samp : processor.defaultdict_accumulator(int) for samp in samples})
        histos['cutflows'] = cutflows
        histos['cutDesc'] = processor.defaultdict_accumulator(str)

        # creating the accumulator
        self._accumulator = processor.dict_accumulator(histos)
        
        # load in cuts module
        self.cutFile = cutFile
        if "/" in self.cutFile: # if cut file is in a different directory
            sys.path.append("/".join(self.cutFile.split("/")[:-1]))
            cutFileName = self.cutFile.split("/")[-1].split(".")[0]
            self.cutLib = importlib.import_module(cutFileName)
            cutList = [c for c in dir(self.cutLib) if "cut" in c]
            cutList = sorted(cutList,key=lambda x: int(x[3:])) # make sure cuts are ordered as they are in the file
            self.cuts = [getattr(self.cutLib,c) for c in cutList]
        else: # cut file is in the same directory (e.g. running on condor)
            cutFileName = self.cutFile.split(".")[0]
            self.cutLib = importlib.import_module(cutFileName)
            cutList = [c for c in dir(self.cutLib) if "cut" in c]
            cutList = sorted(cutList,key=lambda x: int(x[3:])) # make sure cuts are ordered as they are in the file
            self.cuts = [getattr(self.cutLib,c) for c in cutList]


    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self,events):
        samp = events.metadata["dataset"]
        histos = self.accumulator.identity()
        info = self.sampleInfo[samp]

        # Initial number of events
        histos['cutflows'][samp]['initial'] += len(events)

        # Preselection
        routines.selectGoodElesAndVertices(events)

        events.__setitem__("nGoodVtx",ak.count(events.good_RRvtx.vxy,axis=1) + 
                                        ak.count(events.good_LRvtx.vxy,axis=1) +
                                        ak.count(events.good_LLvtx.vxy,axis=1))
        events = events[events.nGoodVtx > 0]

        # pre-computing quantities for cuts
        events.__setitem__("JetMETdPhi",deltaPhi(events.PFJet.corrPhi,events.PFMET.correctedPhi))
        routines.selectBestVertex(events)

        # pre-computing quantities for histograms, as specified in the histo config
        for subroutine in self.subroutines:
            getattr(routines,subroutine)(events)
        
        ###############################
        ######## CUTS & HISTOS ########
        ###############################
        for cut in self.cuts:
            events, cutName, cutDesc, savePlots = cut(events,info)

            histos['cutflows'][samp][cutName] += len(events)
            histos['cutDesc'][cutName] += cutDesc + "@"

            # Fill histograms
            if savePlots:
                self.histoFill(events,histos,samp,cutName)
        
        return histos

    def postprocess(self, accumulator):
        # only need one description per cut name -- adds many during parallel execution
        for cutName in list(accumulator['cutDesc'].keys()):
            accumulator['cutDesc'][cutName] = accumulator['cutDesc'][cutName].split("@")[0]
        return accumulator

class fileAnalyzer(Analyzer):
    def __init__(self,fileList,sample,histoList,routines=[]):
        # load in configs
        if ".root" in fileList:
            self.fileList = [fileList]
        elif type(fileList) == list:
            self.fileList = fileList
        elif ".txt" in fileList:
            with open(fileList) as f:
                self.fileList = f.read().splitlines()
        
        with open(histoList) as f:
            self.histoList = json.load(f)
        self.routines = routines

        self.sample_names = [sample] # list of sample names, readable names are generated from data in the fileList json
        self.sample_locs = {sample : self.fileList} # dictionary mapping sample name to file/directory location
        self.sample_info = {} # dictionary with sample metadata
        self.histos = {}
        self.histoData = {}
        self.nEvents = {}
        self.nEventsProcessed = {}
        self.totalEvents = 0
        
def deltaPhi(v1,v2):
    # copy of the ROOT RVec DeltaPhi function
    # see here https://root.cern/doc/master/RVec_8hxx_source.html#l02742
    M_PI = 3.14159265358979323846264338328
    dPhi = np.fmod(v1-v2,2*M_PI)
    under = ak.values_astype(dPhi < -1*M_PI,np.float32)
    over = ak.values_astype(dPhi > M_PI,np.float32)
    fine = ak.values_astype((dPhi <= M_PI) & (dPhi >= -1*M_PI),np.float32)
    output = fine*dPhi + under*(dPhi + 2.0*M_PI) + over*(dPhi - 2.0*M_PI)
    return output

def loadSchema(fileLoc):
    loc = uproot.open(fileLoc)
    tree = NanoEventsFactory.from_root(loc,treepath="ntuples/outT",schemaclass=MySchema).events()
    return tree

def flatten_fillNone(arr,val):
    return ak.fill_none(ak.flatten(arr),val)

def loadHistoFiles(location):
    files = [f for f in os.listdir(location) if ".coffea" in f]
    histos = {}
    for f in files:
        htemp = coffea.util.load(location+"/"+f)[0]
        for hname in list(htemp.keys()):
            if hname not in histos.keys():
                histos[hname] = htemp[hname].copy()
            else:
                histos[hname] += htemp[hname]
    return histos

def getSampleInfo(histos,hname="ele_kinematics"):
    samps = [s.name for s in histos[hname].axis("sample").identifiers()]
    info = {}
    for s in samps:
        ct = re.findall("ctau-(\d+)",s)[0]
        m, dm = re.findall("Mchi-(\d+p\d)_dMchi-(\d+p\d)",s)[0]
        m = m.replace("p",".")
        dm = dm.replace("p",".")
        entry = "{0}-{1}".format(m,dm)
        if entry not in info.keys():
            info[entry] = {}
        info[entry][ct] = s
    return info

def adjustAxes(ax,hstyle,binName,ct):
    if binName in hstyle.keys():
        style = hstyle[binName]
        if "xrange" in style.keys():
            if ("vxy"  == binName) or ("vz" == binName):
                ax.set_xlim(left=vxy_range[ct][0],right=vxy_range[ct][1])
            elif style['xrange']:
                ax.set_xlim(style['xrange'][0],right=style['xrange'][1])

def rebinHist(hstyle,hist,binName,ct):
    if binName in hstyle.keys():
        style = hstyle[binName]
        if "rebin" in style.keys():
            if ("vxy" == binName) or ("vz" == binName):
                hist = hist.rebin(binName,vxy_rebin[ct])
            if style['rebin']:
                hist = hist.rebin(binName,style['rebin'])
    return hist

def setLogY(ax,hstyle,binName):
    if binName in hstyle.keys():
        style = hstyle[binName]
        if 'log' in style.keys():
            if style['log']:
                ax.set_yscale('log')

def sampToTitle(samp):
    vals = samp.split("_")
    mchi = float(vals[0].split("-")[1])
    dmchi = float(vals[1].split("-")[1])
    ct = float(vals[2].split("-")[1])
    title = r"$m_1 = {0}$ GeV, $m_2 = {1}$ GeV, $c\tau = {2}$ mm".format(int(mchi - dmchi/2),int(mchi+dmchi/2),int(ct))
    return title

def retrieveParams(samp):
    vals = samp.split("_")
    mchi = float(vals[0].split("-")[1])
    dmchi = float(vals[1].split("-")[1])
    ct = float(vals[2].split("-")[1])
    return mchi, dmchi, ct

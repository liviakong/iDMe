from __future__ import with_statement
import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from mySchema import MySchema
from coffea import processor

#from coffea.dataset_tools import (
#    apply_to_fileset,
#    max_chunks,
#    preprocess,
#)
#import dask

import uproot
import awkward as ak
import vector
vector.register_awkward()
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
from collections import defaultdict
from hist import Hist
from hist.axis import StrCategory, Regular, Integer, IntCategory
import hist

match_names = {"Default":"match0","lowpt":"match1"}
vxy_range = {1:[0,20],10:[0,50],100:[0,50],1000:[0,50]}
vxy_rebin = {1:5,10:20,100:20,1000:20}

class Analyzer:
    def __init__(self,fileList,histoList,cuts,max_samples=-1,max_files_per_samp=-1,newCoffea=False):
        # flag to see if we're using new coffea
        self.newCoffea = newCoffea

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
                if 'mZD' in sample['name']:
                    name = "sig_{0}_Mchi-{1}_dMchi-{2}_ctau-{3}_{4}".format(sample['year'],str(sample['Mchi']).replace(".","p"),str(sample['dMchi']).replace(".","p"),sample['ctau'],sample['name'].split('_')[-1])
                else:
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
                if self.newCoffea:
                    self.sample_locs[name] = {'files':{sample['location']:'ntuples/outT'}}
                else:
                    self.sample_locs[name] = [sample['location']]
            elif 'fileset' in sample.keys():
                if self.newCoffea:
                    self.sample_locs[name] = {'files':{f:'ntuples/outT' for f in sample['fileset'] if f.split("/")[-1] not in sample['blacklist']}}
                    if self.max_files_per_samp > 0 and len(self.sample_locs[name]['files']) > self.max_files_per_samp:
                        self.sample_locs[name]['files'] = {k:self.sample_locs[name]['files'][k] for k in list(self.sample_locs[name]['files'].keys())[:self.max_files_per_samp]}
                else:
                    self.sample_locs[name] = [f for f in sample['fileset'] if f.split("/")[-1] not in sample['blacklist']]
                    if self.max_files_per_samp > 0 and len(self.sample_locs[name]) > self.max_files_per_samp:
                        self.sample_locs[name] = self.sample_locs[name][:self.max_files_per_samp]
            else:
                # if the location is a directory, use the xrootd client to get a list of files
                xrdClient = client.FileSystem("root://cmseos.fnal.gov")
                if type(loc) != list:
                    status, flist = xrdClient.dirlist(loc)
                    fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))]
                else:
                    fullList = []
                    for l in loc:
                        status, flist = xrdClient.dirlist(l)
                        fullList.extend(["root://cmsxrootd.fnal.gov/"+l+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))])
                if self.max_files_per_samp > 0 and len(fullList) > self.max_files_per_samp:
                    fullList = fullList[:self.max_files_per_samp]
                if self.newCoffea:
                    self.sample_locs[name] = {'files':{f:'ntuples/outT' for f in fullList}}
                else:
                    self.sample_locs[name] = fullList
            
            self.sample_info[name] = sample
            self.sample_names.append(name)
            loaded += 1

    def process(self,treename='ntuples/outT',execr="iterative",workers=4,dask_client=None,procType='default',**kwargs):
        fileset = self.sample_locs
        if procType == 'default':
            proc = iDMeProcessor(self.sample_names,self.sample_info,self.sample_locs,self.histoFile,self.cuts,mode=self.mode,**kwargs)
        elif procType == 'gen':
            proc = genProcessor(self.sample_names,self.sample_info,self.sample_locs,self.histoFile,self.cuts,mode=self.mode,**kwargs)
        elif procType == 'trig':
            proc = trigProcessor(self.sample_names,self.sample_info,self.sample_locs,self.histoFile,self.cuts,mode=self.mode,**kwargs)
        
        if not self.newCoffea:
            if execr == "iterative":
                executor = processor.IterativeExecutor()
            elif execr == "futures":
                executor = processor.FuturesExecutor(workers=workers)
            elif execr == "dask":
                if dask_client is None:
                    print("Need to supply a dask client!")
                    return
                else:
                    executor = processor.DaskExecutor(client=dask_client)
            else:
                print("Invalid executor type specification!")
                return
            runner = processor.Runner(executor=executor,schema=MySchema,savemetrics=True)
            accumulator = runner(fileset,
                                treename=treename,
                                processor_instance=proc)
        else:
            print("Preprocessing")
            dataset_runnable, dataset_updated = preprocess(fileset,step_size=100_000,files_per_batch=1)
            print("Done Preprocessing")
            to_compute = apply_to_fileset(proc,dataset_runnable,schemaclass=MySchema)
            (accumulator,) = dask.compute(to_compute)
        
        return accumulator

class iDMeProcessor(processor.ProcessorABC):
    def __init__(self,samples,sampleInfo,fileSet,histoFile,cutFile,mode='signal',**kwargs):
        self.samples = samples
        self.sampleInfo = sampleInfo
        self.sampleLocs = fileSet
        self.mode = mode
        
        # load in histogram config
        self.histoMod = importlib.import_module(histoFile)
        self.histoFill = self.histoMod.fillHistos
        self.subroutines = self.histoMod.subroutines
        
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

        self.extraStuff = {}
        for k,v in kwargs.items():
            self.extraStuff[k] = v
            print(f"Registering extra input {k} = {v}")
    
    def process(self,events):
        samp = events.metadata["dataset"]
        info = self.sampleInfo[samp]
        isMC = info["type"] == "signal" or info["type"] == "bkg"
        info['defineGoodVertices'] = routines.defineGoodVertices
        info['selectBestVertex'] = routines.selectBestVertex
        for k,v in self.extraStuff.items():
            info[f"extras_{k}"] = v
        
        #histos = self.histoMod.make_histograms()
        #histos['cutDesc'] = defaultdict(str)
        histObj = self.histoMod.make_histograms(info)
        cutDesc = defaultdict(str)

        cutflow = defaultdict(float)               # efficiency
        cutflow_counts = defaultdict(float)        # xsec-weighted event counts
        cutflow_nevts = defaultdict(int)           # raw event counts
        cutflow_vtx_matched = defaultdict(float)   # (for signal MC) fraction that the selected vertex is truth-matched
        # (for signal MC) also check the above, but only counting the events where both gen ee are reconstructed: dR(reco,gen) < 0.1
        
        if isMC:
            sum_wgt = info["sum_wgt"]
            lumi, unc = getLumi(info['year'])
            xsec = info['xsec']
            # register event weight branch
            events.__setitem__("eventWgt",xsec*lumi*events.genWgt)
        else:
            sum_wgt = info["num_events"]

        # Initial number of events
        if isMC:
            cutflow['all'] += ak.sum(events.genWgt)/sum_wgt
        else:
            cutflow['all'] += len(events)/sum_wgt
        cutflow_nevts['all'] += len(events)

        if info['type'] == "signal":
            cutflow_vtx_matched['all'] += 1 # dummy value before selecting a vertex

            isReconstructed = events.GenEle.matched & events.GenPos.matched
            has_gen_matched_reco_ee_events = events[isReconstructed]
            cutflow_genEEreconstructed['all'] += ak.sum(has_gen_matched_reco_ee_events.genWgt)/sum_wgt
            cutflow_vtx_matched_genEEreconstructed['all'] += 1 # dummy value before selecting a vertex
            
        cutDesc['all'] = 'No cuts@'

        #################################
        ## Calculating Additional Vars ##
        #################################
        events['Electron','mindRj'] = ak.fill_none(ak.min(events.Electron.dRJets,axis=-1),999)
        events['Electron','mindPhiJ'] = ak.fill_none(ak.min(np.abs(events.Electron.dPhiJets),axis=-1),999)
        events['LptElectron','mindRj'] = ak.fill_none(ak.min(events.LptElectron.dRJets,axis=-1),999)
        events['LptElectron','mindPhiJ'] = ak.fill_none(ak.min(np.abs(events.LptElectron.dPhiJets),axis=-1),999)
        events['vtx','mindRj'] = ak.fill_none(ak.min(events.vtx.dRJets,axis=-1),999)
        events['vtx','mindPhiJ'] = ak.fill_none(ak.min(np.abs(events.vtx.dPhiJets),axis=-1),999)
        events['vtx','cosThetaColl'] = (events.vtx.px/events.vtx.pt)*(events.vtx.vx/events.vtx.vxy) + (events.vtx.py/events.vtx.pt)*(events.vtx.vy/events.vtx.vxy)
        events['vtx','vxyCosThetaColl'] = events.vtx.vxy*events.vtx.cosThetaColl
        routines.electronID(events,info) # electron kinematic/ID definition
        routines.jetBtag(events,info['year'])
        if info['type'] == "signal":
            events['GenJetMETdPhi'] = np.abs(deltaPhi(events.GenJet.phi[:,0],events.GenMET.phi))
            events['GenEle','dr'] = events.genEE.dr
            events['GenPos','dr'] = events.genEE.dr
            if "vxy" not in events.genEE.fields:
                events['genEE','vxy'] = events.GenEle.vxy
            routines.genElectronKinematicBins(events)
            #routines.getLptMatchInfoForReg(events)
            routines.genMatchRecoQuantities(events)
        # associate electrons to vertices after all electron-related stuff has been computed
        routines.vtxElectronConnection(events) # associate electrons to vertices
        events['vtx','min_dxy'] = np.minimum(np.abs(events.vtx.e1.dxy),np.abs(events.vtx.e2.dxy))
        events['vtx','eleDphi'] = np.abs(deltaPhi(events.vtx.e1.phi,events.vtx.e2.phi))
        #if info['type'] == 'signal':
            #routines.genMatchExtraVtxVariables(events)
        
        #################################
        ##### Hard-coded basic cuts #####
        #################################
        # 1 or 2 jets in the event
        nJets = ak.count(events.PFJet.pt,axis=1)
        events = events[(nJets>0) & (nJets<3)]
        # needs a good vertex
        routines.defineGoodVertices(events,version='none') # define "good" vertices based on whether associated electrons pass ID cuts
        events = events[events.nGoodVtx > 0]
        # define "selected" vertex based on selection criteria in the routine (nominally: lowest chi2)
        routines.selectBestVertex(events)
        #events = routines.selectTrueVertex(events,events.good_vtx)

        # Fill cutflow after baseline selection
        if isMC:
            cutflow['hasVtx'] += ak.sum(events.genWgt)/sum_wgt
        else:
            cutflow['hasVtx'] += len(events)/sum_wgt
        cutflow_nevts['hasVtx'] += len(events)
        cutDesc['hasVtx'] = 'Baseline Selection@'

        # For signal, (1) check if the vertex ee are gen-matched (2) check if the event has ee that are gen-matched
        if info['type'] == "signal":
            vtx_matched_events = events[events.sel_vtx.isMatched]
            cutflow_vtx_matched['hasVtx'] += ak.sum(vtx_matched_events.genWgt)/ak.sum(events.genWgt)
            
            isReconstructed = events.GenEle.matched & events.GenPos.matched
            has_gen_matched_reco_ee_events = events[isReconstructed]
            cutflow_genEEreconstructed['hasVtx'] += ak.sum(has_gen_matched_reco_ee_events.genWgt)/sum_wgt
            vtx_matched_events_genEEreconstructed = has_gen_matched_reco_ee_events[has_gen_matched_reco_ee_events.sel_vtx.isMatched]
            cutflow_vtx_matched_genEEreconstructed['hasVtx'] += ak.sum(vtx_matched_events_genEEreconstructed.genWgt)/ak.sum(has_gen_matched_reco_ee_events.genWgt)
        
        # Compute miscellaneous extra variables -- add anything you want to this function
        # don't need this for now - can activate if need be
        #routines.miscExtraVariables(events)
        #if info['type'] == "signal":
        #    routines.miscExtraVariablesSignal(events)

        # computing any extra quantities specified in the histogram config file
        for subroutine in self.subroutines:
            getattr(routines,subroutine)(events)
        
        ###############################
        ######## CUTS & HISTOS ########
        ###############################
        for cut in self.cuts:
            events, cutName, cutDescription, savePlots = cut(events,info)
            if isMC:
                cutflow[cutName] += ak.sum(events.genWgt)/sum_wgt
            else:
                cutflow[cutName] += len(events)/sum_wgt
            cutflow_nevts[cutName] += len(events)            
            if info['type'] == "signal":
                vtx_matched_events = events[events.sel_vtx.isMatched]
                cutflow_vtx_matched[cutName] += ak.sum(vtx_matched_events.genWgt)/ak.sum(events.genWgt)

                isReconstructed = events.GenEle.matched & events.GenPos.matched
                has_gen_matched_reco_ee_events = events[isReconstructed]
                cutflow_genEEreconstructed[cutName] += ak.sum(has_gen_matched_reco_ee_events.genWgt)/sum_wgt
                vtx_matched_events_genEEreconstructed = has_gen_matched_reco_ee_events[has_gen_matched_reco_ee_events.sel_vtx.isMatched]
                cutflow_vtx_matched_genEEreconstructed[cutName] += ak.sum(vtx_matched_events_genEEreconstructed.genWgt)/ak.sum(has_gen_matched_reco_ee_events.genWgt)
            cutDesc[cutName] += cutDescription + "@"

            # Fill histograms
            if savePlots and len(events) > 0: # fixes some bugginess trying to fill histograms with empty arrays
                self.histoFill(events,histObj,samp,cutName,info,sum_wgt=sum_wgt)
        
        for k in cutflow.keys():
            if isMC:
                cutflow_counts[k] = xsec*lumi*cutflow[k]
                cutflow_counts_genEEreconstructed[k] = xsec*lumi*cutflow_genEEreconstructed[k]
            else:
                cutflow_counts[k] = sum_wgt*cutflow[k]
                cutflow_counts_genEEreconstructed = sum_wgt*cutflow_genEEreconstructed[k]
        
        histos = histObj.histograms
        histos['cutDesc'] = cutDesc
        histos['cutflow'] = {samp:cutflow}
        histos['cutflow_cts'] = {samp:cutflow_counts}
        histos['cutflow_nevts'] = {samp:cutflow_nevts}
        histos['cutflow_vtx_matched'] = {samp:cutflow_vtx_matched}

        return histos

    def postprocess(self, accumulator):
        # only need one description per cut name -- adds many during parallel execution
        for cutName in list(accumulator['cutDesc'].keys()):
            accumulator['cutDesc'][cutName] = accumulator['cutDesc'][cutName].split("@")[0]
        return accumulator

class genProcessor(iDMeProcessor):
    def process(self,events):
        samp = events.metadata["dataset"]
        info = self.sampleInfo[samp]
        
        #histos = self.histoMod.make_histograms()
        #histos['cutDesc'] = defaultdict(str)
        histObj = self.histoMod.make_histograms(info)
        
        cutDesc = defaultdict(str)
        cutflow = defaultdict(float)               # efficiency
        cutflow_counts = defaultdict(float)        # xsec-weighted event counts
        cutflow_nevts = defaultdict(int)           # raw event counts
        
        sum_wgt = info["sum_wgt"]
        lumi, unc = getLumi(info['year'])
        xsec = info['xsec']

        # register event weight branch
        events.__setitem__("eventWgt",xsec*lumi*events.genWgt)

        # Initial number of events
        cutflow['all'] += ak.sum(events.genWgt)/sum_wgt
        cutflow_nevts['all'] += len(events)
        cutDesc['all'] = 'No cuts'

        #################################
        ## Calculating Additional Vars ##
        #################################
        events['Electron','mindRj'] = ak.fill_none(ak.min(events.Electron.dRJets,axis=-1),999)
        events['Electron','mindPhiJ'] = ak.fill_none(ak.min(events.Electron.dPhiJets,axis=-1),999)
        events['LptElectron','mindRj'] = ak.fill_none(ak.min(events.LptElectron.dRJets,axis=-1),999)
        events['LptElectron','mindPhiJ'] = ak.fill_none(ak.min(events.LptElectron.dPhiJets,axis=-1),999)
        events['vtx','mindRj'] = ak.fill_none(ak.min(events.vtx.dRJets,axis=-1),999)
        events['vtx','mindPhiJ'] = ak.fill_none(ak.min(events.vtx.dPhiJets,axis=-1),999)
        routines.electronID(events,info) # electron kinematic/ID definition
        if info['type'] == "signal":
            events['GenJetMETdPhi'] = deltaPhi(events.GenJet.phi[:,0],events.GenMET.phi)
            events['GenEle','dr'] = events.genEE.dr
            events['GenPos','dr'] = events.genEE.dr
            if "vxy" not in events.genEE.fields:
                events['genEE','vxy'] = events.GenEle.vxy
            routines.genElectronKinematicBins(events)
            routines.getLptMatchInfoForReg(events)
            routines.genMatchRecoQuantities(events)
        # associate electrons to vertices after all electron-related stuff has been computed
        routines.vtxElectronConnection(events) # associate electrons to vertices
        routines.defineGoodVertices(events) # define "good" vertices based on whether associated electrons pass ID cuts
        if info['type'] == 'signal':
            routines.genMatchExtraVtxVariables(events)

        # initial histogram fill
        self.histoFill(events,histObj,samp,"no_presel",info,sum_wgt=sum_wgt)

        #################################
        #### Hard-coded basic cuts ######
        #################################
        # 1 or 2 jets in the event
        nJets = ak.count(events.PFJet.pt,axis=1)
        events = events[(nJets>0) & (nJets<3)]
        
        #################################
        #### Demand >= 1 ee vertices ####
        #################################
        events.__setitem__("nGoodVtx",ak.count(events.good_vtx.vxy,axis=1))
        events = events[events.nGoodVtx > 0]
        # define "selected" vertex based on selection criteria in the routine (nominally: lowest chi2)
        routines.selectBestVertex(events)

        # Fill cutflow after baseline selection
        cutflow['hasVtx'] += ak.sum(events.genWgt)/sum_wgt
        cutflow_nevts['hasVtx'] += len(events)
        cutDesc['hasVtx'] = 'Baseline Selection'
        
        # Compute miscellaneous extra variables -- add anything you want to this function
        routines.miscExtraVariables(events)

        # computing any extra quantities specified in the histogram config file
        for subroutine in self.subroutines:
            getattr(routines,subroutine)(events)
        
        ###############################
        ######## CUTS & HISTOS ########
        ###############################
        for cut in self.cuts:
            events, cutName, cutDescription, savePlots = cut(events,info)
            cutflow[cutName] += ak.sum(events.genWgt)/sum_wgt
            cutflow_nevts[cutName] += len(events)            
            cutDesc[cutName] += cutDescription + "@"
            # Fill histograms
            if savePlots and len(events) > 0: # fixes some bugginess trying to fill histograms with empty arrays
                self.histoFill(events,histObj,samp,cutName,info,sum_wgt=sum_wgt)
        
        for k in cutflow.keys():
            cutflow_counts[k] = xsec*lumi*cutflow[k]

        histos = histObj.histograms
        histos['cutDesc'] = cutDesc
        histos['cutflow'] = {samp:cutflow}
        histos['cutflow_cts'] = {samp:cutflow_counts}
        histos['cutflow_nevts'] = {samp:cutflow_nevts}
        
        return histos

class trigProcessor(iDMeProcessor):
    def process(self,events):
        samp = events.metadata["dataset"]

        info = self.sampleInfo[samp]
        isMC = info["type"] == "signal" or info["type"] == "bkg"
        if isMC:
            sum_wgt = info["sum_wgt"]
            lumi, unc = getLumi(info['year'])
            xsec = info['xsec']
            # register event weight branch
            events.__setitem__("eventWgt",xsec*lumi*events.genWgt)
        else:
            sum_wgt = info["num_events"]

        # fill histo
        if isMC:
            events['wgt'] = events.eventWgt/sum_wgt
        else:
            events['wgt'] = 1
        
        # define histo
        MET_passTrig = Hist(StrCategory([],name="samp",label="Sample Name",growth=True),
                               Regular(220,80,300,name="met",label="met"),
                               IntCategory([0,1],name="passTrig",label="passTrig"),
                               storage=hist.storage.Weight())

        MET_all = Hist(StrCategory([],name="samp",label="Sample Name",growth=True),
                               Regular(300,0,300,name="met",label="met"),
                               storage=hist.storage.Weight())

        jet_pt_passTrig = Hist(StrCategory([],name="samp",label="Sample Name",growth=True),
                               Regular(500,0,500,name="pt",label="pt]"),
                               IntCategory([0,1],name="passTrig",label="passTrig"),
                               storage=hist.storage.Weight())

        jet_pt_all = Hist(StrCategory([],name="samp",label="Sample Name",growth=True),
                               Regular(500,0,500,name="pt",label="pt]"),
                               storage=hist.storage.Weight())

        # full MET histo
        MET_all.fill(samp=samp,
                            met=events.PFMET.pt,
                            weight=events.wgt)

        # require other stuff
        nJets = ak.count(events.PFJet.pt,axis=1)
        events = events[nJets>0]
        # fill jet pt histo
        jet_pt_all.fill(samp=samp,pt=events.PFJet.pt[:,0],weight=events.wgt)

        # iDM-like selection
        lead_pt = events.PFJet.pt[:,0]
        lead_eta = np.abs(events.PFJet.eta[:,0])
        ele_cut = (events.Electron.pt > 30) & (events.Electron.IDcutLoose==1)
        lpt_ele_cut = (events.LptElectron.pt > 30)
        n_ele_cut = ak.any(ele_cut,axis=1) | ak.any(lpt_ele_cut,axis=1)
        cut = (lead_pt > 80) & (lead_eta < 2.4) & (n_ele_cut)
        events = events[cut]
        
        # require reference trigger
        #events = events[events.trig.HLT_PFJet80==1]

        # fill histos
        MET_passTrig.fill(samp=samp,
                             met=events.PFMET.pt,
                             passTrig=ak.values_astype(events.trig.HLT_PFMET120_PFMHT120_IDTight,int),
                             weight=events.wgt)
        jet_pt_passTrig.fill(samp=samp,pt=events.PFJet.pt[:,0],passTrig=ak.values_astype(events.trig.HLT_PFMET120_PFMHT120_IDTight,int),weight=events.wgt)

        output = {"MET_passTrig":MET_passTrig,"MET_all":MET_all, "jet_pt_all":jet_pt_all, "jet_pt_passTrig":jet_pt_passTrig}

        return output

    def postprocess(self, accumulator):
        return accumulator

class fileSkimmer:
    def __init__(self,sampFile,sampleInfo,cutFile,mode='signal'):
        self.sampleInfo = sampleInfo
        self.sampFile = sampFile
        fname = sampFile.split("/")[-1]
        self.outFileName = fname.replace(".root","_skimmed.root")
        self.mode = mode
        
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
    
    def skim(self):
        with uproot.open(self.sampFile) as input_file:
            events = NanoEventsFactory.from_root(input_file,treepath="ntuples/outT",schemaclass=MySchema).events()
            info = self.sampleInfo
            sum_wgt = info["sum_wgt"]
            lumi, unc = getLumi(info['year'])
            xsec = self.sampleInfo['xsec']

            # register event weight branch
            events.__setitem__("eventWgt",xsec*lumi*events.genWgt/sum_wgt)
            # Preselection
            routines.selectGoodElesAndVertices(events)
            events.__setitem__("nGoodVtx",ak.count(events.good_vtx.vxy,axis=1))
            events = events[events.nGoodVtx > 0]
            
            # pre-computing quantities for cuts
            routines.selectBestVertex(events)
            # computing some signal-only diagnostic quantities
            if info['type'] == "signal":
                e1_match = routines.matchedVertexElectron(events,1)
                e2_match = routines.matchedVertexElectron(events,2)
                events["sel_vtx","match"] = ak.values_astype(ak.where(e1_match*e2_match == -1,2,ak.where(np.abs(e1_match)+np.abs(e2_match) > 0,1,0)),np.int32)
            else:
                events["sel_vtx","match"] = ak.zeros_like(events.sel_vtx.pt,dtype=np.int32)

            ###############################
            ######## CUTS & HISTOS ########
            ###############################
            for cut in self.cuts:
                events, cutName, cutDesc, savePlots = cut(events,info)
            
            if len(events.genWgt) == 0:
                pass
            else:
                output_tree = {}
                vtx_vars = [f for f in events.sel_vtx.fields if f!="e1" and f!="e2"]
                for v in vtx_vars:
                    if v == "e1_typ" or v == "e2_typ":
                        output_tree[f"sel_vtx_{v}"] = ak.Array(ak.where(events.sel_vtx[v]=="R",1,2).to_list())
                    elif v == "typ":
                        output_tree[f"sel_vtx_{v}"] = ak.Array(ak.where(events.sel_vtx[v]=="RR",1,ak.where(events.sel_vtx[v]=="LR",2,3)).to_list())
                    else:
                        output_tree[f"sel_vtx_{v}"] = ak.Array(events.sel_vtx[v].to_list())
                ele_fields = events.sel_vtx.e1.fields
                for ef in ele_fields:
                    output_tree[f"sel_e1_{ef}"] = ak.Array(events.sel_vtx.e1[ef].to_list())
                    output_tree[f"sel_e2_{ef}"] = ak.Array(events.sel_vtx.e2[ef].to_list())
                output_tree['genWgt'] = ak.Array(events.genWgt.to_list())
                output_tree['eventWgt'] = ak.Array(events.eventWgt.to_list())
                additional_fields = {
                                        'CaloMET':['ET','pt','phi'],
                                        'Photon':['et','eta','phi'],
                                        'PFMET':['ET','pt','phi'],
                                        'PFJet':['pt','eta','phi','bTag','METdPhi'],
                                        'Electron':["*"],
                                        'LptElectron':["*"],
                                        'vtx':["*"]
                                    }
                for af in additional_fields.keys():
                    subfields = additional_fields[af]
                    if subfields == ["*"]:
                        subfields = events[af].fields
                    allvars = {}
                    for subf in subfields:
                        if af == "vtx" and ("e1" in subf or "e2" in subf):
                            continue
                        if af == "vtx" and subf=="typ":
                            allvars[subf] = ak.Array(ak.where(events[af][subf]=="RR",1,ak.where(events[af][subf]=="LR",2,3)).to_list())
                        else:
                            allvars[subf] = ak.Array(events[af][subf].to_list())
                    output_tree[af] = ak.zip(allvars)
                with uproot.recreate(self.outFileName) as outfile:
                    outfile['outT'] = output_tree

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

def getLumi(year):
    # recommendations from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
    lumi, unc = 0, 0
    if year == 2016:
        lumi = 36.31
        unc = 0.012*lumi # 1.2 percent
    if year == 2017:
        lumi = 41.48
        unc = 0.023*lumi # 2.3 percent
    if year == 2018:
        lumi = 59.83
        unc = 0.025*lumi # 2.5 percent
    return lumi, unc

def loadSchema(fileLoc):
    loc = uproot.open(fileLoc)
    tree = NanoEventsFactory.from_root(loc,treepath="ntuples/outT",schemaclass=MySchema).events()
    return tree

def loadNano(fileLoc):
    loc = uproot.open(fileLoc)
    tree = NanoEventsFactory.from_root(loc,treepath="ntuples/outT",schemaclass=NanoAODSchema).events()
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

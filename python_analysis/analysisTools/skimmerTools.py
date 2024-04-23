from __future__ import with_statement
import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from mySchema import MySchema
from coffea import processor
from coffea.nanoevents.methods import vector
from coffea.processor import column_accumulator
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
from collections import defaultdict
from analysisTools import getLumi, deltaPhi
import plotTools as ptools

class Skimmer:
    def __init__(self,fileList,cuts,max_samples=-1,max_files_per_samp=-1):
        # load in file config
        if type(fileList) == str and ".json" in fileList:
            with open(fileList) as f:
                self.fileList = json.load(f)
        else:
            self.fileList = fileList

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
                self.sample_locs[name] = [f for f in sample['fileset'] if f.split("/")[-1] not in sample['blacklist']]
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
                if self.max_files_per_samp > 0:
                    fullList = fullList[:self.max_files_per_samp] if len(fullList) > self.max_files_per_samp else fullList
                self.sample_locs[name] = fullList
            
            self.sample_info[name] = sample
            self.sample_names.append(name)
            loaded += 1

    def process(self,treename='ntuples/outT',execr="iterative",workers=4):
        fileset = self.sample_locs
        proc = makeBDTInputs(self.sample_names,self.sample_info,self.sample_locs,self.cuts,mode=self.mode)
        if execr == "iterative":
            executor = processor.IterativeExecutor()
        elif execr == "futures":
            executor = processor.FuturesExecutor(workers=workers)
        else:
            print("Invalid executor type specification!")
            exit
        runner = processor.Runner(executor=executor,schema=MySchema,savemetrics=True)
        accumulator = runner(fileset,
                            treename=treename,
                            processor_instance=proc)
        return accumulator

class makeBDTInputs(processor.ProcessorABC):
    def __init__(self,samples,sampleInfo,fileSet,cutFile,mode='signal'):
        self.samples = samples
        self.sampleInfo = sampleInfo
        self.sampleLocs = fileSet
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
    
    def process(self,events):
        samp = events.metadata["dataset"]
        outputs = {}
        info = self.sampleInfo[samp]
        sum_wgt = info["sum_wgt"]
        lumi, unc = getLumi(info['year'])
        xsec = info['xsec']

        # register event weight branch
        events.__setitem__("eventWgt",xsec*lumi*events.genWgt)
        events.__setitem__("eventWgtNorm",xsec*lumi*events.genWgt/sum_wgt)

        if info['type'] == "signal":
            events.__setitem__("m1",ptools.signalPoint(samp)['m1'])
            events.__setitem__("delta",ptools.signalPoint(samp)['delta'])
            events.__setitem__("ctau",ptools.signalPoint(samp)['ctau'])

        #################################
        #### Hard-coded basic cuts ######
        #################################
        # 1 or 2 jets in the event
        nJets = ak.count(events.PFJet.pt,axis=1)
        events = events[(nJets>0) & (nJets<3)]

        #################################
        ## Calculating Additional Vars ##
        #################################
        routines.electronJetSeparation(events) # dR and dPhi between electrons and jets
        routines.electronIsoConePtSum(events) # for each electron, compute pT sum of any other electrons in event within dR < 0.3 of it
        routines.electronID(events) # electron kinematic/ID definition
        routines.vtxElectronConnection(events) # associate electrons to vertices
        routines.defineGoodVertices(events) # define "good" vertices based on whether associated electrons pass ID cuts
        
        #################################
        #### Demand >= 1 ee vertices ####
        #################################
        events.__setitem__("nGoodVtx",ak.count(events.good_vtx.vxy,axis=1))
        events = events[events.nGoodVtx > 0]
        # define "selected" vertex based on selection criteria in the routine (nominally: lowest chi2)
        #routines.selectBestVertex(events)

        if info['type'] == "signal":
            events = routines.selectTrueVertex(events, events.good_vtx)
        else:
            routines.selectBestVertex(events)

        # Compute miscellaneous extra variables -- add anything you want to this function
        routines.miscExtraVariables(events)
        if info['type'] == "signal":
            routines.miscExtraVariablesSignal(events)

        ###############################
        ######## CUTS & HISTOS ########
        ###############################
        for cut in self.cuts:
            events, cutName, cutDesc, savePlots = cut(events,info) # apply cuts on jets, BDT will be only used for vertex related variables
        
        # filling outputs dictionary
        e1 = events.sel_vtx.e1
        e2 = events.sel_vtx.e2

        deltadxy = np.abs(np.abs(e1.dxy) - np.abs(e2.dxy))
        mindxy = np.minimum(np.abs(e1.dxy),np.abs(e2.dxy))
        maxdxy = np.maximum(np.abs(e1.dxy),np.abs(e2.dxy))
        meandxy = np.mean(np.abs(e1.dxy),np.abs(e2.dxy))
        abs_deta = np.abs(e1.eta - e2.eta)
        abs_dphi = np.abs(e1.phi - e2.phi)
        
        outputs['wgt'] = column_accumulator(events["eventWgt"].to_numpy())
        outputs['wgt_norm'] = column_accumulator(events["eventWgtNorm"].to_numpy())
        outputs['lead_jet_eta'] = column_accumulator(events.PFJet.eta[:,0].to_numpy())
        outputs['lead_jet_pt'] = column_accumulator(events.PFJet.pt[:,0].to_numpy())
        outputs['jetMETdPhi'] = column_accumulator(np.abs(events.PFJet.METdPhi[:,0]).to_numpy())
        outputs['minJetMETdPhi'] = column_accumulator(ak.min(np.abs(events.PFJet.METdPhi),axis=1).to_numpy())
        outputs['sel_vtx_sign'] = column_accumulator(events.sel_vtx.sign.to_numpy())
        outputs['sel_vtx_chi2'] = column_accumulator(events.sel_vtx.reduced_chi2.to_numpy())
        outputs['sel_vtx_METdPhi'] = column_accumulator(np.abs(events.sel_vtx.METdPhi).to_numpy())
        outputs['sel_vtx_m'] = column_accumulator(events.sel_vtx.m.to_numpy())
        outputs['sel_vtx_dR'] = column_accumulator(events.sel_vtx.dR.to_numpy())
        outputs['sel_vtx_minDxy'] = column_accumulator(np.minimum(np.abs(events.sel_vtx.e1.dxy),np.abs(events.sel_vtx.e2.dxy)).to_numpy())
        outputs['met_leadPt_ratio'] = column_accumulator((events.PFMET.pt/events.PFJet.pt[:,0]).to_numpy())
        outputs["vxy"] = column_accumulator(events.sel_vtx.vxy.to_numpy())
        outputs["vxy_signif"] = column_accumulator((events.sel_vtx.vxy/events.sel_vtx.sigmavxy).to_numpy())
        outputs["sel_e1_dxy"] = column_accumulator(np.abs(e1.dxy).to_numpy())
        outputs["sel_e1_dxySignif"] = column_accumulator((np.abs(e1.dxy)/e1.dxyErr).to_numpy())
        outputs["sel_e2_dxy"] = column_accumulator(np.abs(e2.dxy).to_numpy())
        outputs["sel_e2_dxySignif"] = column_accumulator((np.abs(e2.dxy)/e2.dxyErr).to_numpy())
        if info['type'] == 'signal':
            outputs['sel_vtx_match'] = column_accumulator(events.sel_vtx.match.to_numpy())
            outputs['m1'] = column_accumulator(events["m1"].to_numpy())
            outputs['delta'] = column_accumulator(events["delta"].to_numpy())
            outputs['ctau'] = column_accumulator(events["ctau"].to_numpy())

        outputs["sel_vtx_cos_collinear"] = column_accumulator(events.cos_collinear.to_numpy())
        outputs["sel_vtx_prod_eta"] = column_accumulator((e1.eta*e2.eta).to_numpy())
        outputs["sel_vtx_sign_prod_eta"] = column_accumulator(((e1.eta*e2.eta)/np.abs(e1.eta*e2.eta)).to_numpy())

        outputs["sel_vtx_pt"] = column_accumulator(events.sel_vtx.pt.to_numpy())
        outputs["sel_vtx_projectedLxy"] = column_accumulator(events.projectedLxy.to_numpy())
        outputs["sel_vtx_pt_e1_over_pt_e2"] = column_accumulator((np.minimum(e1.pt, e2.pt)/np.maximum(e1.pt, e2.pt)).to_numpy())
        outputs["sel_vtx_pt_over_m"] = column_accumulator((events.sel_vtx.pt/events.sel_vtx.m).to_numpy())

        outputs["delta_dxy_over_mindxy"] = column_accumulator((deltadxy/mindxy).to_numpy())
        outputs["delta_dxy_over_maxdxy"] = column_accumulator((deltadxy/maxdxy).to_numpy())
        outputs["delta_dxy_over_meandxy"] = column_accumulator((deltadxy/meandxy).to_numpy())

        outputs["delta_eta_over_delta_phi"] = column_accumulator((abs_deta/abs_dphi).to_numpy())
        outputs["log_delta_eta_over_delta_phi"] = column_accumulator(np.log10(abs_deta/abs_dphi).to_numpy())
        
        return outputs

    def postprocess(self, accumulator):
        return accumulator
    
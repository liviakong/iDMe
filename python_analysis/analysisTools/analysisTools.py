from __future__ import with_statement
import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
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
NanoAODSchema.warn_missing_crossrefs = False

match_names = {"Default":"match0","lowpt":"match1"}
vxy_range = {1:[0,20],10:[0,50],100:[0,50],1000:[0,50]}
vxy_rebin = {1:5,10:20,100:20,1000:20}

class Analyzer:
    def __init__(self,fileList,histoList,histoStyle=None):
        # load in configs
        with open(fileList) as f:
            self.fileList = json.load(f)
        with open(histoList) as f:
            self.histoList = json.load(f)
        if histoStyle != None:
            with open(histoStyle) as f:
                self.histoStyle = json.load(f)

        self.sample_names = [] # list of sample names, readable names are generated from data in the fileList json
        self.sample_locs = {} # dictionary mapping sample name to file/directory location
        self.sample_info = {} # dictionary with sample metadata
        self.histos = {}
        self.histoData = {}
        self.histoStyle = {}
        self.nEvents = {}
        self.nEventsProcessed = {}
        self.totalEvents = 0
        
        self.loadFiles()
        self.loadHistos()
    
    def loadFiles(self):
        for sample in self.fileList:
            mode = sample['type']
            if mode == 'signal':
                name = "sig_{0}_Mchi-{1:.1f}_dMchi-{2:.1f}_ctau-{3}".format(sample['year'],sample['Mchi'],sample['dMchi'],sample['ctau'])
            elif mode == 'bkg':
                name = "bkg_{0}_{1}".format(sample['year'],sample['name'])
            elif mode == 'data':
                name = "data_{0}_{1}".format(sample['year'],sample['name'])
            loc = sample['location']
            if '.root' in loc:
                # if the location is just a single file, load it in
                self.sample_locs[name] = [sample['location']]
            else:
                # if the location is a directory, use the xrootd client to get a list of files
                xrdClient = client.FileSystem("root://cmseos.fnal.gov")
                status, flist = xrdClient.dirlist(loc)
                fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if '.root' in item.name]
                self.sample_locs[name] = fullList
            self.sample_info[name] = sample
            self.sample_names.append(name)

    def loadTable(self,fileLoc,treeName,treeDir='ntuples_gbm/'):
        NanoAODSchema.warn_missing_crossrefs = False
        loc = uproot.open(fileLoc)
        tree = NanoEventsFactory.from_root(loc,treepath=treeDir+treeName,schemaclass=NanoAODSchema).events()
        nEntries = tree.num_entries
        return tree, nEntries
    
    def loadHistos(self):
        for histo in self.histoList:
            cats = histo['categories']
            bins = histo['bins']
            vname = histo['varname']
            self.histoData[vname] = {'bins':[],'cats':[]}
            axes = []
            for k in cats.keys():
                axes.append(hist.Cat(k,cats[k]))
                self.histoData[vname]['cats'].append(k)
            for b in bins:
                if (len(b) == 3): # third argument is a list of bins
                    axes.append(hist.Bin(b[0],b[1],b[2]))
                if (len(b) == 5):
                    axes.append(hist.Bin(b[0],b[1],b[2],b[3],b[4]))
                self.histoData[vname]['bins'].append(b[0])
            h = hist.Hist("Events",axes=axes)
            self.histos[vname] = h

    def getHistoStyle(self,key):
        if key in self.histoStyle.keys():
            return self.histoStyle[key]
        else:
            return {}

class iDMeProcessor(processor.ProcessorABC):
    def __init__(self,samples,sampleInfo,fileSet,histos,mode='signal'):
        self.samples = samples
        self.sampleInfo = sampleInfo
        self.sampleLocs = fileSet
        self.mode = mode
        self._accumulator = processor.dict_accumulator(histos)

    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self,events):
        samp = events.metadata["dataset"]
        histos = self.accumulator.identity()
        eles = events.Electron

        ele_kin = histos['ele_kinematics']
        ele_trkHits = histos['ele_trackHits']
        ele_trkQual = histos['ele_trackQual']

        match_ele_kin = histos['matched_ele_kinematics']
        match_ele_gen_kin = histos['matched_ele_gen_kinematics']
        match_ele_trkHits = histos['matched_ele_trackHits']
        match_ele_trkQual = histos['matched_ele_trackQual']
        match_ele_gen_disp = histos['matched_ele_gen_displacement']

        match_mtx = histos['match_matrix']
        conv_match_mtx = histos['conv_match_matrix']

        gen_ele_disp = histos['gen_displacement']
        
        regreg_vtx_vxy = histos["regreg_reco_vertex_vxy"]
        regreg_vtx_stats = histos["regreg_reco_vertex_stats"]
        lowlow_vtx_vxy = histos["lowlow_reco_vertex_vxy"]
        lowlow_vtx_stats = histos["lowlow_reco_vertex_stats"]
        lowreg_vtx_vxy = histos["lowreg_reco_vertex_vxy"]
        lowreg_vtx_stats = histos["lowreg_reco_vertex_stats"]
        
        ################ FILLING START ################    
        # Filling basic electron kinematics
        eles = events.Electron
        lpt_eles = events.LptElectron
        convs = events.Conversion
        eles_dir = ak.zip(
            {
                "x":eles.eta,
                "y":eles.phi
            },
            with_name="TwoVector"
        )
        lpt_eles_dir = ak.zip(
            {
                "x":lpt_eles.eta,
                "y":lpt_eles.phi
            },
            with_name="TwoVector"
        )
        conv1_dir = ak.zip(
            {
                "x":convs.trk1_innerEta,
                "y":convs.trk1_innerPhi,
            },
            with_name="TwoVector"
        )
        conv2_dir = ak.zip(
            {
                "x":convs.trk2_innerEta,
                "y":convs.trk2_innerPhi,
            },
            with_name="TwoVector"
        )

        all_eles_dir = ak.concatenate([eles_dir,lpt_eles_dir,conv1_dir,conv2_dir],axis=1)
        all_eles_labels = ak.concatenate([ak.ones_like(eles.eta),2*ak.ones_like(lpt_eles.eta),
                                          3*ak.ones_like(convs.trk1_innerEta),3*ak.ones_like(convs.trk2_innerEta)],axis=1)

        ele_kin.fill(sample=samp,
                     ele_type="Default",
                     pt=ak.flatten(eles.pt),
                     eta=ak.flatten(eles.eta),
                     phi=ak.flatten(eles.phi))
        ele_kin.fill(sample=samp,
                     ele_type="lowpt",
                     pt=ak.flatten(lpt_eles.pt),
                     eta=ak.flatten(lpt_eles.eta),
                     phi=ak.flatten(lpt_eles.phi))
        
        ele_trkHits.fill(sample=samp,
                         ele_type="Default",
                         numTrkHits=ak.flatten(eles.numTrackerHits),
                         numPixHits=ak.flatten(eles.numPixHits),
                         numStripHits=ak.flatten(eles.numStripHits))
        ele_trkHits.fill(sample=samp,
                         ele_type="lowpt",
                         numTrkHits=ak.flatten(lpt_eles.numTrackerHits),
                         numPixHits=ak.flatten(lpt_eles.numPixHits),
                         numStripHits=ak.flatten(lpt_eles.numStripHits))

        ele_trkQual.fill(sample=samp,
                         ele_type="Default",
                         chi2=ak.flatten(eles.trkChi2),
                         trkIso=ak.flatten(eles.trkIso))
        ele_trkQual.fill(sample=samp,
                         ele_type="lowpt",
                         chi2=ak.flatten(lpt_eles.trkChi2),
                         trkIso=ak.flatten(lpt_eles.trkIso))
        
        ################ Gen-matching ################
        if self.mode == "signal":
            # Match 0 - default electrons matched to gen
            gen_eles = events.GenPart[np.abs(events.GenPart.ID) == 11]
            gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])
            gen_pos = ak.flatten(events.GenPart[events.GenPart.ID == -11])
            gen_ele_dir = ak.zip(
                {
                    "x":gen_ele.eta,
                    "y":gen_ele.phi,
                },
                with_name="TwoVector"
            )
            gen_pos_dir = ak.zip(
                {
                    "x":gen_pos.eta,
                    "y":gen_pos.phi,
                },
                with_name="TwoVector"
            )

            dR_e = (all_eles_dir - gen_ele_dir).r
            dR_p = (all_eles_dir - gen_pos_dir).r
            nearest_e = ak.argmin(dR_e,axis=1,keepdims=True)
            nearest_p = ak.argmin(dR_p,axis=1,keepdims=True)
            mindR_e = ak.fill_none(ak.flatten(dR_e[nearest_e]),999)
            mindR_p = ak.fill_none(ak.flatten(dR_p[nearest_p]),999)
            type_e_nearest = ak.fill_none(ak.flatten(all_eles_labels[nearest_e]),0)
            type_p_nearest = ak.fill_none(ak.flatten(all_eles_labels[nearest_p]),0)
            
            match_e = mindR_e < 0.1
            match_p = mindR_p < 0.1
            type_e_nearest = type_e_nearest*ak.values_astype(match_e,int) # sets match type to 0 (no match) where dR < 0.1 criterion not satisfied
            type_p_nearest = type_p_nearest*ak.values_astype(match_p,int)
            match_mtx.fill(sample=samp,
                           e_matchType=type_e_nearest,
                           p_matchType=type_p_nearest)

            ## Conversion-specific studies
            d1e = (conv1_dir - gen_ele_dir).r
            d1p = (conv1_dir - gen_pos_dir).r
            d2e = (conv2_dir - gen_ele_dir).r
            d2p = (conv2_dir - gen_pos_dir).r
            d1e = d1e[ak.argmin(d1e,axis=1,keepdims=True)]
            d1p = d1p[ak.argmin(d1p,axis=1,keepdims=True)]
            d2e = d2e[ak.argmin(d2e,axis=1,keepdims=True)]
            d2p = d2p[ak.argmin(d2p,axis=1,keepdims=True)]
            match_1e2p_in = (d1e < 0.1) & (d2p < 0.1) & (d1p > 0.1) & (d2e > 0.1)
            match_1p2e_in = (d1p < 0.1) & (d2e < 0.1) & (d1e > 0.1) & (d2p > 0.1)
            fullUniqueMatch_in = match_1e2p_in | match_1p2e_in
            fullUniqueMatch_in = ak.fill_none(ak.flatten(fullUniqueMatch_in),False)
            fullMatch_in = ((d1e < 0.1) | (d1p < 0.1)) & ((d2e < 0.1) | (d2p < 0.1))
            fullMatch_in = ak.fill_none(ak.flatten(fullMatch_in),False)
            any_match_in = (d1e < 0.1) | (d1p < 0.1) | (d2e < 0.1) | (d2p < 0.1)
            any_match_in = ak.fill_none(ak.flatten(any_match_in),False)

            ele_kin.fill(sample=samp,ele_type="Generator",
                        pt=ak.flatten(gen_eles.pt),eta=ak.flatten(gen_eles.eta),phi=ak.flatten(gen_eles.phi))
            
            d1_match0 = gen_eles[:,0].delta_r(eles)
            d2_match0 = gen_eles[:,1].delta_r(eles)
            match0 = np.logical_or(d1_match0 < 0.1,d2_match0 < 0.1)
            match0_gen_indices = ak.values_astype(d1_match0 > d2_match0,int) # select index 1 (True) if d2 < d1 (gen ele at index 1 is nearer)
            eles_match0 = eles[match0]
            eles_match0_gen = gen_eles[match0_gen_indices][match0]
            match_ele_kin.fill(sample=samp,
                               match="Default",
                               pt=ak.flatten(eles_match0.pt),
                               eta=ak.flatten(eles_match0.eta),
                               phi=ak.flatten(eles_match0.phi))
            match_ele_gen_kin.fill(sample=samp,
                                   match="Default",
                                   pt=ak.flatten(eles_match0_gen.pt),
                                   eta=ak.flatten(eles_match0_gen.eta),
                                   phi=ak.flatten(eles_match0_gen.phi))
            match_ele_trkHits.fill(sample=samp,
                                   match="Default",
                                   numTrkHits=ak.flatten(eles_match0.numTrackerHits),
                                   numPixHits=ak.flatten(eles_match0.numPixHits),
                                   numStripHits=ak.flatten(eles_match0.numStripHits),
                                   vxy=ak.flatten(eles_match0_gen.vxy))
            match_ele_trkQual.fill(sample=samp,
                                   match="Default",
                                   chi2=ak.flatten(eles_match0.trkChi2),
                                   trkIso=ak.flatten(eles_match0.trkIso))
            match_ele_gen_disp.fill(sample=samp,
                                    match="Default",
                                    vxy=ak.flatten(eles_match0_gen.vxy),
                                    vz=ak.flatten(eles_match0_gen.vz))
            
            # Match 1 - low-pT electrons matched to gen
            d1_match1 = gen_eles[:,0].delta_r(lpt_eles)
            d2_match1 = gen_eles[:,1].delta_r(lpt_eles)
            match1 = np.logical_or(d1_match1 < 0.1,d2_match1 < 0.1)
            match1_gen_indices = ak.values_astype(d1_match1 > d2_match1,int) # select index 1 (True) if d2 < d1 (gen ele at index 1 is nearer)
            eles_match1 = lpt_eles[match1]
            eles_match1_gen = gen_eles[match1_gen_indices][match1]
            match_ele_kin.fill(sample=samp,
                               match="lowpt",
                               pt=ak.flatten(eles_match1.pt),
                               eta=ak.flatten(eles_match1.eta),
                               phi=ak.flatten(eles_match1.phi))
            match_ele_gen_kin.fill(sample=samp,
                                   match="lowpt",
                                   pt=ak.flatten(eles_match1_gen.pt),
                                   eta=ak.flatten(eles_match1_gen.eta),
                                   phi=ak.flatten(eles_match1_gen.phi))
            match_ele_trkHits.fill(sample=samp,
                                   match="lowpt",
                                   numTrkHits=ak.flatten(eles_match1.numTrackerHits),
                                   numPixHits=ak.flatten(eles_match1.numPixHits),
                                   numStripHits=ak.flatten(eles_match1.numStripHits),
                                   vxy=ak.flatten(eles_match1_gen.vxy))
            match_ele_trkQual.fill(sample=samp,
                                   match="lowpt",
                                   chi2=ak.flatten(eles_match1.trkChi2),
                                   trkIso=ak.flatten(eles_match1.trkIso))
            match_ele_gen_disp.fill(sample=samp,
                                    match="lowpt",
                                    vxy=ak.flatten(eles_match1_gen.vxy),
                                    vz=ak.flatten(eles_match1_gen.vz))

        ################ Reco Electron Vertex Histograms ################
        vtx_regreg = events.EleVertexRegReg
        vtx_lowlow = events.EleVertexLowLow
        vtx_lowreg = events.EleVertexLowReg
        # Selecting only the entries with a valid vertex, i.e. vxy > 0
        vtx_regreg = vtx_regreg[vtx_regreg.vxy > 0]
        vtx_lowlow = vtx_lowlow[vtx_lowlow.vxy > 0]
        vtx_lowreg = vtx_lowreg[vtx_lowreg.vxy > 0]
        # Filters for selecting opposite charge dielectron vertices
        opp_regreg = vtx_regreg.sign == -1
        opp_lowlow = vtx_lowlow.sign == -1
        opp_lowreg = vtx_lowreg.sign == -1
        # reg-reg, i.e. two regular electrons matched to a vertex
        regreg_vtx_vxy.fill(sample=samp,
                            sign="all",
                            vxy=ak.flatten(vtx_regreg.vxy),
                            sigma_vxy=ak.flatten(vtx_regreg.sigmavxy))
        regreg_vtx_vxy.fill(sample=samp,
                            sign="opp",
                            vxy=ak.flatten(vtx_regreg[opp_regreg].vxy),
                            sigma_vxy=ak.flatten(vtx_regreg[opp_regreg].sigmavxy))
        regreg_vtx_stats.fill(sample=samp,
                              sign="all",
                              chi2=ak.flatten(vtx_regreg.reduced_chi2),
                              signif=ak.flatten(vtx_regreg.vxy/vtx_regreg.sigmavxy))
        regreg_vtx_stats.fill(sample=samp,
                              sign="opp",
                              chi2=ak.flatten(vtx_regreg[opp_regreg].reduced_chi2),
                              signif=ak.flatten(vtx_regreg[opp_regreg].vxy/vtx_regreg[opp_regreg].sigmavxy))
        # low-low, i.e. two low-pT electrons matched to a vertex
        lowlow_vtx_vxy.fill(sample=samp,
                            sign="all",
                            vxy=ak.flatten(vtx_lowlow.vxy),
                            sigma_vxy=ak.flatten(vtx_lowlow.sigmavxy))
        lowlow_vtx_vxy.fill(sample=samp,
                            sign="opp",
                            vxy=ak.flatten(vtx_lowlow[opp_lowlow].vxy),
                            sigma_vxy=ak.flatten(vtx_lowlow[opp_lowlow].sigmavxy))
        lowlow_vtx_stats.fill(sample=samp,
                              sign="all",
                              chi2=ak.flatten(vtx_lowlow.reduced_chi2),
                              signif=ak.flatten(vtx_lowlow.vxy/vtx_lowlow.sigmavxy))
        lowlow_vtx_stats.fill(sample=samp,
                              sign="opp",
                              chi2=ak.flatten(vtx_lowlow[opp_lowlow].reduced_chi2),
                              signif=ak.flatten(vtx_lowlow[opp_lowlow].vxy/vtx_lowlow[opp_lowlow].sigmavxy))
        # low-reg, i.e. one low-pT and one default electron matched to a vertex
        lowreg_vtx_vxy.fill(sample=samp,
                            sign="all",
                            vxy=ak.flatten(vtx_lowreg.vxy),
                            sigma_vxy=ak.flatten(vtx_lowreg.sigmavxy))
        lowreg_vtx_vxy.fill(sample=samp,
                            sign="opp",
                            vxy=ak.flatten(vtx_lowreg[opp_lowreg].vxy),
                            sigma_vxy=ak.flatten(vtx_lowreg[opp_lowreg].sigmavxy))
        lowreg_vtx_stats.fill(sample=samp,
                              sign="all",
                              chi2=ak.flatten(vtx_lowreg.reduced_chi2),
                              signif=ak.flatten(vtx_lowreg.vxy/vtx_lowreg.sigmavxy))
        lowreg_vtx_stats.fill(sample=samp,
                              sign="opp",
                              chi2=ak.flatten(vtx_lowreg[opp_lowreg].reduced_chi2),
                              signif=ak.flatten(vtx_lowreg[opp_lowreg].vxy/vtx_lowreg[opp_lowreg].sigmavxy))
            
        ################ GenParticle Histograms ################
        gen_ele_disp.fill(sample=samp,
                          vxy=ak.flatten(gen_eles.vxy),
                          vz=ak.flatten(gen_eles.vz))
        
        return histos

    def postprocess(self, accumulator):
        return accumulator

def plotHistos(alzr,density=False):
    plt.style.use('plt_settings.mplstyle')
    for hname in alzr.histos.keys():
        h = alzr.histos[hname]
        hstyle = alzr.histoStyle[hname]
        axes = h.axes()
        cats = {}
        bins = []
        for a in axes:
            if type(a) == coffea.hist.hist_tools.Cat:
                cats[a.name] = [idt.name for idt in a.identifiers()]
            if type(a) == coffea.hist.hist_tools.Bin:
                bins.append(a.name)
        if 'sample' not in cats.keys():
            print("Histogram {0} has no 'sample' category, skipping".format(hname))
        for samp in cats['sample']:
            # load sample histogram and parameters 
            hsamp = h.integrate("sample",int_range=samp)
            mchi,dmchi,ct = retrieveParams(samp)
            if density:
                plot_dir = "plots/Mchi-{0}_dMchi-{1}/normalized/".format(mchi,dmchi)
            else:
                plot_dir = "plots/Mchi-{0}_dMchi-{1}/plots/".format(mchi,dmchi)
            if not os.path.isdir(plot_dir):
                os.mkdir(plot_dir)
            title = sampToTitle(samp)
            overlay = None
            title_prefix = ""
            if "ele_type" in cats.keys():
                overlay = "ele_type"
                title_prefix = "compare_eles_"
            if "match" in cats.keys():
                overlay = "match"
                title_prefix = "compare_match_"
            
            # loop over 1D variables
            for b in bins:
                # for each variable, integrate out the others
                hsampbin = hsamp
                others = [p for p in bins if p != b]
                for o in others:
                    hsampbin = hsampbin.integrate(o)
                
                #make the plot
                plt.figure(figsize=(10,8))
                #custom styling
                if hstyle[b]["rebin"] : hsampbin = hsampbin.rebin(b,hstyle[b]["rebin"])
                if ("vxy" in b) or ("vz" in b):
                    hsampbin = hsampbin.rebin(b,vxy_rebin[ct])
                if hstyle[b]["log"] : plt.yscale("log")
                hist.plot1d(hsampbin,overlay=overlay,density=density)
                if hstyle[b]["xrange"] : plt.xlim(hstyle[b]["xrange"])
                if ("vxy" in b) or ("vz" in b):
                    plt.xlim(vxy_range[ct])
                plt.title(title)
                plt.savefig(plot_dir+title_prefix+b+"_ct-"+str(int(ct))+".pdf")
                plt.close()

def plotEfficiencies(alzr):
    plt.style.use('plt_settings.mplstyle')
    # Gen-matching efficiencies
    for samp in alzr.samples.keys():
        mchi,dmchi,ct = retrieveParams(samp)
        title = sampToTitle(samp)
        plot_dir = "plots/Mchi-{0}_dMchi-{1}/efficiencies/".format(mchi,dmchi)
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)
        # Gen-matched electron efficiencies, kinematics
        gen_kinematics = alzr.histos["ele_kinematics"].integrate("sample",int_range=samp).integrate("ele_type","Generator")
        matched_gen_kinematics = alzr.histos["matched_ele_gen_kinematics"].integrate("sample",int_range=samp)
        kin_bins = alzr.histoData["ele_kinematics"]['bins']
        
        gen_disp = alzr.histos["gen_displacement"].integrate("sample",int_range=samp)
        matched_gen_disp = alzr.histos["matched_ele_gen_displacement"].integrate("sample",int_range=samp)
        disp_bins = alzr.histoData["gen_displacement"]['bins']
        
        hstyle_kin = alzr.histoStyle["ele_kinematics"]
        hstyle_disp = alzr.histoStyle["gen_displacement"]

        match_cats = [a.name for a in matched_gen_kinematics.axis("match").identifiers()]
        for cat in match_cats:                
            for b in kin_bins:
                hnum_kin = matched_gen_kinematics.integrate("match",int_range=cat)
                hden_kin = gen_kinematics
                other = [bn for bn in kin_bins if bn != b]
                for b2 in other:
                    hnum_kin = hnum_kin.integrate(b2)
                    hden_kin = hden_kin.integrate(b2)
                plt.figure(figsize=(10,8))
                #custom styling
                if hstyle_kin[b]["rebin"] : 
                    hnum_kin = hnum_kin.rebin(b,hstyle_kin[b]["rebin"])
                    hden_kin = hden_kin.rebin(b,hstyle_kin[b]["rebin"])
                if hstyle_kin[b]["log"] : plt.yscale("log")
                hist.plotratio(hnum_kin,hden_kin,ax=plt.gca(),error_opts={'color': 'k', 'marker': '.'},unc='poisson-ratio')
                if hstyle_kin[b]["xrange"] : plt.xlim(hstyle_kin[b]["xrange"])
                match_cat = "Match Category : {0}".format(cat)
                plt.title(title)
                plt.text(0,0.95,match_cat,horizontalalignment='left',transform=plt.gca().transAxes)
                plt.savefig(plot_dir+match_names[cat]+"_"+b+"_eff_ct-"+str(int(ct))+".pdf")
                plt.close()

            for b in disp_bins:
                hnum_disp = matched_gen_disp.integrate("match",int_range=cat)
                hden_disp = gen_disp
                other = [bn for bn in disp_bins if bn != b]
                for b2 in other:
                    hnum_disp = hnum_disp.integrate(b2)
                    hden_disp = hden_disp.integrate(b2)
                plt.figure(figsize=(10,8))
                #custom styling
                if hstyle_disp[b]["log"] : plt.yscale("log")
                hnum_disp = hnum_disp.rebin(b,vxy_rebin[ct])
                hden_disp = hden_disp.rebin(b,vxy_rebin[ct])
                hist.plotratio(hnum_disp,hden_disp,ax=plt.gca(),error_opts={'color': 'k', 'marker': '.'},unc='poisson-ratio')
                plt.xlim(vxy_range[ct])
                match_cat = "Match Category : {0}".format(cat)
                plt.title(title)
                plt.text(0,0.95,match_cat,horizontalalignment='left',transform=plt.gca().transAxes)
                plt.savefig(plot_dir+match_names[cat]+"_"+b+"_eff_ct-"+str(int(ct))+".pdf")
                plt.close()

def compare1DHists(an1,an2,label1,label2,file_prefix="",plot_base=""):
    if plot_base == "":
        print("Please specify a directory to save the outputs! Use plot_base=[plot_dir] when calling the function")
        return
    
    subdirs = plot_base.split("/")
    for k in range(len(subdirs)):
        d = "/".join(subdirs[:k+1])
        if not os.path.isdir(d):
            os.mkdir(d)
        
    hlist = list(an1.histos.keys())
    samples = list(an1.samples.keys())

    for samp in samples:
        m,dm,ct = retrieveParams(samp)
        title = sampToTitle(samp)
        dirname = "Mchi-{0}_dMchi-{1}_ctau-{2}".format(m,dm,ct)
        if not os.path.isdir(plot_base+dirname): 
            os.mkdir(plot_base+dirname)
        
        # Loop through histograms
        for hname in hlist:
            pref = plot_base+dirname+"/"
            if not os.path.isdir(pref+hname):
                os.mkdir(pref+hname)
            pref = pref + hname + "/"
            
            h1 = an1.histos[hname].integrate("sample",samp)
            h2 = an2.histos[hname].integrate("sample",samp)
            hStyle = an1.getHistoStyle(hname)
            hData = an1.histoData[hname]
            categories = hData['cats']
            categories = [c for c in categories if c != 'sample']
            bins = hData['bins']
            subcategories = {c:[b.name for b in h1.axis(c).identifiers()] for c in categories}
            
            # Loop through numerical bins (e.g. pT, eta, phi, etc.)
            for b in bins:
                # Loop through categories (e.g. match category, etc.)
                for cat in categories:
                    other_bins = [bn for bn in bins if bn != b]
                    other_cats = [c for c in categories if c != cat]
                    subcats = subcategories[cat]
                    # Loop through 'bins' (sub-categories) of the selected category (e.g. match0, match1, etc.)
                    for subcat in subcats:
                        # Pick out sub-category
                        h1_plot = h1.integrate(cat,subcat)
                        h2_plot = h2.integrate(cat,subcat)
                        # Integrate over other categories and axes
                        h1_plot = h1_plot.sum(*(other_cats+other_bins),overflow='all')
                        h2_plot = h2_plot.sum(*(other_cats+other_bins),overflow='all')
                        # Plot the histograms
                        plt.style.use('plt_settings.mplstyle')
                        fig,ax = plt.subplots(figsize=(10,8))
                        h1_plot = rebinHist(hStyle,h1_plot,b,ct)
                        h2_plot = rebinHist(hStyle,h2_plot,b,ct)
                        setLogY(ax,hStyle,b)
                        hist.plot1d(h1_plot)
                        hist.plot1d(h2_plot)
                        adjustAxes(ax,hStyle,b,ct)
                        ax.legend([label1,label2])
                        plt.title(title)
                        out_fname = pref+file_prefix+b+'_'+cat+'_'+subcat+'.pdf'
                        plt.savefig(out_fname)
                        plt.close()

def compare1DHists_difference(an1,an2,file_prefix="",plot_base=""):
    if plot_base == "":
        print("Please specify a directory to save the outputs! Use plot_base=[plot_dir] when calling the function")
        return
    
    subdirs = plot_base.split("/")
    for k in range(len(subdirs)):
        d = "/".join(subdirs[:k+1])
        if not os.path.isdir(d):
            os.mkdir(d)
        
    hlist = list(an1.histos.keys())
    samples = list(an1.samples.keys())

    for samp in samples:
        m,dm,ct = retrieveParams(samp)
        title = sampToTitle(samp)
        dirname = "Mchi-{0}_dMchi-{1}_ctau-{2}".format(m,dm,ct)
        if not os.path.isdir(plot_base+dirname): 
            os.mkdir(plot_base+dirname)
        
        # Loop through histograms
        for hname in hlist:
            pref = plot_base+dirname+"/"
            if not os.path.isdir(pref+hname):
                os.mkdir(pref+hname)
            pref = pref + hname + "/"
            
            h1 = an1.histos[hname].integrate("sample",samp)
            h2 = an2.histos[hname].integrate("sample",samp)
            hStyle = an1.getHistoStyle(hname)
            hData = an1.histoData[hname]
            categories = hData['cats']
            categories = [c for c in categories if c != 'sample']
            bins = hData['bins']
            subcategories = {c:[b.name for b in h1.axis(c).identifiers()] for c in categories}
            
            # Loop through numerical bins (e.g. pT, eta, phi, etc.)
            for b in bins:
                # Loop through categories (e.g. match category, etc.)
                for cat in categories:
                    other_bins = [bn for bn in bins if bn != b]
                    other_cats = [c for c in categories if c != cat]
                    subcats = subcategories[cat]
                    # Loop through 'bins' (sub-categories) of the selected category (e.g. match0, match1, etc.)
                    for subcat in subcats:
                        # Pick out sub-category
                        h1_plot = h1.integrate(cat,subcat)
                        h2_plot = h2.integrate(cat,subcat)
                        # Integrate over other categories and axes
                        h1_plot = h1_plot.sum(*(other_cats+other_bins),overflow='all')
                        h2_plot = h2_plot.sum(*(other_cats+other_bins),overflow='all')
                        # Plot the histograms
                        plt.style.use('plt_settings.mplstyle')
                        fig,ax = plt.subplots(figsize=(10,8))
                        h1_plot = rebinHist(hStyle,h1_plot,b,ct)
                        h2_plot = rebinHist(hStyle,h2_plot,b,ct)
                        h2_plot.scale(-1)
                        h1_plot.add(h2_plot)
                        hist.plot1d(h1_plot,line_opts={})
                        plt.axhline(y=0,color='gray',linestyle='--')
                        ax.get_legend().remove()
                        adjustAxes(ax,hStyle,b,ct)
                        counts = list(h1_plot.values().values())[0]
                        plt.ylim([1.5*np.min(counts),1.5*np.max(counts)])
                        plt.title(title)
                        out_fname = pref+file_prefix+b+'_'+cat+'_'+subcat+'.pdf'
                        plt.savefig(out_fname)
                        plt.close()

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
            
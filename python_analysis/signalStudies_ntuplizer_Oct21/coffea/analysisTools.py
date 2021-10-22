import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
import coffea.hist as hist
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import json
import os

match_names = {"Default":"match0",r"Low $p_T$":"match1"}
vxy_range = {1:[0,2],10:[0,20],100:[0,50],1000:[0,50]}
vxy_rebin = {1:1,10:10,100:20,1000:20}

class Analyzer:
    def __init__(self,fileList,histoList,histoStyle):
        self.samples = {}
        self.histos = {}
        self.files = {}
        self.recoTs = {}
        self.genTs = {}
        self.histoData = {}
        
        with open(fileList) as f:
            js = json.load(f)
        for sample in js:
            path = sample['file']
            name = "Mchi-{0:.1f}_dMchi-{1:.1f}_ctau-{2}".format(sample['Mchi'],sample['dMchi'],sample['ctau'])
            self.samples[name] = path
            
        with open(histoList) as f:
            js = json.load(f)
        for histo in js:
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

        with open(histoStyle) as f:
            self.histoStyle = json.load(f)
            
    def loadTables(self):
        for s in self.samples.keys():
            file = uproot.open(self.samples[s])
            self.files[s] = file
            self.recoTs[s] = NanoEventsFactory.from_root(file,treepath='ntuples_gbm/recoT',schemaclass=NanoAODSchema).events()
            self.genTs[s] = NanoEventsFactory.from_root(file,treepath='ntuples_gbm/genT',schemaclass=NanoAODSchema).events()
            
            
    def fillHistos(self):
        ele_kin = self.histos['ele_kinematics']
        ele_trkHits = self.histos['ele_trackHits']
        ele_trkQual = self.histos['ele_trackQual']

        match_ele_kin = self.histos['matched_ele_kinematics']
        match_ele_gen_kin = self.histos['matched_ele_gen_kinematics']
        match_ele_trkHits = self.histos['matched_ele_trackHits']
        match_ele_trkQual = self.histos['matched_ele_trackQual']
        match_ele_gen_disp = self.histos['matched_ele_gen_displacement']

        gen_ele_disp = self.histos['gen_displacement']
        
        
        for samp in self.samples.keys():
            recoT = self.recoTs[samp]
            genT = self.genTs[samp]
            
            # Filling basic electron kinematics
            eles = recoT.Electron
            lpt_eles = recoT.LptElectron
            gen_eles = genT.GenPart[np.abs(genT.GenPart.ID) == 11]

            ele_kin.fill(sample=samp,ele_type="Default",
                         pt=ak.flatten(eles.pt),eta=ak.flatten(eles.eta),phi=ak.flatten(eles.phi))
            ele_kin.fill(sample=samp,ele_type=r"Low $p_T$",
                         pt=ak.flatten(lpt_eles.pt),eta=ak.flatten(lpt_eles.eta),phi=ak.flatten(lpt_eles.phi))
            ele_kin.fill(sample=samp,ele_type="Generator",
                         pt=ak.flatten(gen_eles.pt),eta=ak.flatten(gen_eles.eta),phi=ak.flatten(gen_eles.phi))
            
            ele_trkHits.fill(sample=samp,ele_type="Default",
                              numTrkHits=ak.flatten(eles.numTrackerHits),numPixHits=ak.flatten(eles.numPixHits),numStripHits=ak.flatten(eles.numStripHits))
            ele_trkHits.fill(sample=samp,ele_type=r"Low $p_T$",
                              numTrkHits=ak.flatten(lpt_eles.numTrackerHits),numPixHits=ak.flatten(lpt_eles.numPixHits),numStripHits=ak.flatten(lpt_eles.numStripHits))

            ele_trkQual.fill(sample=samp,ele_type="Default",
                              chi2=ak.flatten(eles.trkChi2),trkIso=ak.flatten(eles.trkIso))
            ele_trkQual.fill(sample=samp,ele_type=r"Low $p_T$",
                              chi2=ak.flatten(lpt_eles.trkChi2),trkIso=ak.flatten(lpt_eles.trkIso))
            
            ################ Gen-matching ################
            # Match 0 - default electrons matched to gen
            d1_match0 = gen_eles[:,0].delta_r(eles)
            d2_match0 = gen_eles[:,1].delta_r(eles)
            match0 = np.logical_or(d1_match0 < 0.1,d2_match0 < 0.1)
            match0_gen_indices = ak.values_astype(d1_match0 > d2_match0,int) # select index 1 (True) if d2 < d1 (gen ele at index 1 is nearer)
            eles_match0 = eles[match0]
            eles_match0_gen = gen_eles[match0_gen_indices][match0]
            match_ele_kin.fill(sample=samp,match="Default",
                               pt=ak.flatten(eles_match0.pt),eta=ak.flatten(eles_match0.eta),phi=ak.flatten(eles_match0.phi))
            match_ele_gen_kin.fill(sample=samp,match="Default",
                                   pt=ak.flatten(eles_match0_gen.pt),eta=ak.flatten(eles_match0_gen.eta),phi=ak.flatten(eles_match0_gen.phi))
            match_ele_trkHits.fill(sample=samp,match="Default",
                                    numTrkHits=ak.flatten(eles_match0.numTrackerHits),numPixHits=ak.flatten(eles_match0.numPixHits),numStripHits=ak.flatten(eles_match0.numStripHits))
            match_ele_trkQual.fill(sample=samp,match="Default",
                                    chi2=ak.flatten(eles_match0.trkChi2),trkIso=ak.flatten(eles_match0.trkIso))
            match_ele_gen_disp.fill(sample=samp,match="Default",
                                    vxy=ak.flatten(eles_match0_gen.vxy),vz=ak.flatten(eles_match0_gen.vz))
            
            # Match 1 - low-pT electrons matched to gen
            d1_match1 = gen_eles[:,0].delta_r(lpt_eles)
            d2_match1 = gen_eles[:,1].delta_r(lpt_eles)
            match1 = np.logical_or(d1_match1 < 0.1,d2_match1 < 0.1)
            match1_gen_indices = ak.values_astype(d1_match1 > d2_match1,int) # select index 1 (True) if d2 < d1 (gen ele at index 1 is nearer)
            eles_match1 = lpt_eles[match1]
            eles_match1_gen = gen_eles[match1_gen_indices][match1]
            match_ele_kin.fill(sample=samp,match=r"Low $p_T$",
                               pt=ak.flatten(eles_match1.pt),eta=ak.flatten(eles_match1.eta),phi=ak.flatten(eles_match1.phi))
            match_ele_gen_kin.fill(sample=samp,match=r"Low $p_T$",
                                   pt=ak.flatten(eles_match1_gen.pt),eta=ak.flatten(eles_match1_gen.eta),phi=ak.flatten(eles_match1_gen.phi))
            match_ele_trkHits.fill(sample=samp,match=r"Low $p_T$",
                                    numTrkHits=ak.flatten(eles_match1.numTrackerHits),numPixHits=ak.flatten(eles_match1.numPixHits),numStripHits=ak.flatten(eles_match1.numStripHits))
            match_ele_trkQual.fill(sample=samp,match=r"Low $p_T$",
                                    chi2=ak.flatten(eles_match1.trkChi2),trkIso=ak.flatten(eles_match1.trkIso))
            match_ele_gen_disp.fill(sample=samp,match=r"Low $p_T$",
                                    vxy=ak.flatten(eles_match1_gen.vxy),vz=ak.flatten(eles_match1_gen.vz))

            ################ Reco Electron Vertex Histograms ################
            vtx_regreg = recoT.EleVertex.regreg
            vtx_lowlow = recoT.EleVertex.lowlow
            vtx_lowreg = recoT.EleVertex.lowreg


            ################ GenParticle Histograms ################
            gen_ele_disp.fill(sample=samp,vxy=ak.flatten(gen_eles.vxy),vz=ak.flatten(gen_eles.vz))

    def plotHistos(self,density=False):
        plt.style.use('plt_settings.mplstyle')
        for hname in self.histos.keys():
            h = self.histos[hname]
            hstyle = self.histoStyle[hname]
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

    def plotEfficiencies(self):
        plt.style.use('plt_settings.mplstyle')
        # Gen-matching efficiencies
        for samp in self.samples.keys():
            mchi,dmchi,ct = retrieveParams(samp)
            title = sampToTitle(samp)
            plot_dir = "plots/Mchi-{0}_dMchi-{1}/efficiencies/".format(mchi,dmchi)
            if not os.path.isdir(plot_dir):
                os.mkdir(plot_dir)
            # Gen-matched electron efficiencies, kinematics
            gen_kinematics = self.histos["ele_kinematics"].integrate("sample",int_range=samp).integrate("ele_type","Generator")
            matched_gen_kinematics = self.histos["matched_ele_gen_kinematics"].integrate("sample",int_range=samp)
            kin_bins = self.histoData["ele_kinematics"]['bins']
            
            gen_disp = self.histos["gen_displacement"].integrate("sample",int_range=samp)
            matched_gen_disp = self.histos["matched_ele_gen_displacement"].integrate("sample",int_range=samp)
            disp_bins = self.histoData["gen_displacement"]['bins']
            
            hstyle_kin = self.histoStyle["ele_kinematics"]
            hstyle_disp = self.histoStyle["gen_displacement"]

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
            
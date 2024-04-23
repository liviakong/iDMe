import hist
import coffea.util as util
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re
import os
import sys
import mplhep as hep
import utils

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

def makeCDF(h,start,stop,bins=100,right=True,nevents=False,category=False):
    x = np.linspace(start,stop,bins)
    n_tot = h.sum(flow=True).value
    if not category:
        if right:
            yields = np.array([h[complex(f"{xi}j")::sum].value for xi in x])
        else:
            yields = np.array([h[:complex(f"{xi}j"):sum].value for xi in x])
    else:
        edges_ordered = h.axes[0].edges[::-1] if right else h.axes[0].edges
        yields_raw = np.cumsum(h.values()[::-1]) if right else np.cumsum(h.values())
        # add extra point to yields to draw steps
        yields = []
        x = []
        for i,y in enumerate(yields_raw):
            yields.extend([y,y])
            x.extend([edges_ordered[i],edges_ordered[i+1]])
        yields = np.array(yields)
        x = np.array(x)
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

# Plot efficiency type stuff
def plot_signal_efficiency(sig_histo, df, m1s, deltas, ctaus, doLog = True, ylabel = '', title = ''):
    cuts = utils.get_signal_list_of_cuts(sig_histo)

    m1_list = []
    for point in df.index.values:
        sig_dict = signalPoint(point)
        m1 = int(sig_dict['m1'])
        m1_list.append(m1)

    df['m1'] = m1_list
    df = df.sort_values(by=['m1']) # sort by m1
    df.pop('m1')
        
    for point in df.index.values:
        sig_dict = signalPoint(point)
        m1 = int(sig_dict['m1'])
        delta = sig_dict['delta']
        dmchi = sig_dict['dmchi']
        ctau = int(sig_dict['ctau'])
        
        if (m1 in m1s) and (delta in deltas):
            if ctau in ctaus:
                plt.plot(cuts, df.loc[point], label=f'({m1}, {dmchi}) GeV, ctau = {int(ctau)}mm')

    plt.grid()

    if doLog:
        plt.yscale('log')
    
    plt.ylabel(ylabel)
    plt.title(title)
    
    plt.xticks(ticks = np.arange(len(cuts)), labels = cuts, rotation = 45, ha = 'right')
    
    plt.legend()
    plt.show()

def plot_bkg_efficiency(bkg_histos, df, doLog = True, ylabel = '', title = ''):
    processes = df.index.values.tolist()
    cuts = utils.get_bkg_list_of_cuts(bkg_histos)

    # Color map for each process
    # cmap = mpl.colormaps['Set3'].colors
    # cmap = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"] # cms-recommended for 6-color scheme
    cmap = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"] # cms-recommended
    
    colors = { 'W+jets': cmap[0],
               'Z+jets': cmap[1],
               'QCD': cmap[2],
               'DY': cmap[3],
               'Top': cmap[4],
               'TTJetsDiLept': cmap[5],
               'Diboson': cmap[6],
               'Triboson': cmap[7],
               'Total': cmap[8]
    }

    for process in processes:
        plt.plot(cuts, df.loc[process], label=process, color = colors[process])

    plt.grid()

    if doLog:
        plt.yscale('log')
    
    plt.ylabel(ylabel)
    plt.title(title)
    
    plt.xticks(ticks = np.arange(len(cuts)), labels = cuts, rotation = 45, ha = 'right')
    
    plt.legend()
    plt.show()

# Plot kinematics
def plot_signal_1D(ax, sig_histo, m1, delta, ctau, plot_dict, style_dict):
    # get signal point info
    si = utils.get_signal_point_dict(sig_histo)
    samp_df = si[(si.m1 == m1) & (si.delta == delta) & (si.ctau == ctau)]
    
    samp = samp_df.name[0]

    m1 = samp_df.m1[0]
    dmchi = samp_df.dmchi[0]
    ctau = samp_df.ctau[0]
    label = f'({m1}, {dmchi}) GeV, ctau = {int(ctau)}mm'

    if style_dict['label'] != None:
        label = style_dict['label']
    
    # get histogram from coffea output
    histo = sig_histo[plot_dict['variable']][{"samp":samp, "cut": plot_dict['cut']}]

    # rebinning
    histo = histo[::style_dict['rebin']]

    # set x range manually
    if style_dict['xlim'] != None:
        xlim = style_dict['xlim']
        xbin_range = np.where((histo.axes.edges[0] > xlim[0]) & (histo.axes.edges[0] < xlim[1]))[0]
        histo = histo[ int(xbin_range[0])-1:int(xbin_range[-1]+1) ]

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])
    else:    
        binwidth = histo.axes.widths[0][0]
        if style_dict['doDensity']:
            ax.set_ylabel(f'A.U./{binwidth:.3f}')
        else:
            ax.set_ylabel(f'Events/{binwidth:.3f}')

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')

    # Plot
    hep.histplot(histo, yerr=style_dict['doYerr'], density=style_dict['doDensity'], ax=ax, histtype='step', flow=style_dict['flow'], label = label)


def plot_signal_2D(ax, sig_histo, m1, delta, ctau, plot_dict, style_dict):
    # get signal point info
    si = utils.get_signal_point_dict(sig_histo)
    samp_df = si[(si.m1 == m1) & (si.delta == delta) & (si.ctau == ctau)]
    
    samp = samp_df.name[0]

    m1 = samp_df.m1[0]
    dmchi = samp_df.dmchi[0]
    ctau = samp_df.ctau[0]
    label = f'({m1}, {dmchi}) GeV, ctau = {int(ctau)}mm'
    
    # get histogram from coffea output
    histo = sig_histo[plot_dict['variable']][{"samp":samp, "cut": plot_dict['cut']}]

    # rebinning
    histo = histo[::style_dict['xrebin'],::style_dict['yrebin']]

    # set x range manually
    if style_dict['xlim'] != None:
        xlim = style_dict['xlim']
        xbin_range = np.where((histo.axes.edges[0] > xlim[0]) & (histo.axes.edges[0] < xlim[1]))[0]
        histo = histo[ int(xbin_range[0])-1:int(xbin_range[-1]+1), : ]
    if style_dict['ylim'] != None:
        ylim = style_dict['ylim']
        ybin_range = np.where((histo.axes.edges[1] > ylim[0]) & (histo.axes.edges[1] < ylim[1]))[1]
        histo = histo[ :, int(ybin_range[0]):int(ybin_range[-1]+1) ]
    
    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])
    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')

    if style_dict['doLogz']:
        hep.hist2dplot(histo, flow=style_dict['flow'], norm=mpl.colors.LogNorm(), ax=ax)
    else:
        hep.hist2dplot(histo, flow=style_dict['flow'], ax=ax)


def plot_bkg_1d(ax, bkg_histos, plot_dict, style_dict, processes = 'all'):  

    if processes == 'all':
        #processes = bkg_histos.keys()

        list_cut_index = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=True)
        list_cut_name = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=False)
        
        cut_name = plot_dict['cut']
        
        df = utils.get_bkg_cutflow_df(bkg_histos, 'cutflow_cts').iloc[:-1]
        
        df = df[list_cut_name[list_cut_index.index(cut_name)]]
        
        processes = df.index[df != 0].to_list()
    
    # if process is given as a list, i.e. ['DY', 'W+jets'], plot only these processes in the list; otherwise, plot all as default

    bkg={}
    bkg[plot_dict['variable']] = {process:bkg_histos[process][plot_dict['variable']][{"samp":sum}] for process in processes}
    
    # sort the histograms by the entries and stack
    for process in processes:
        entries = {process: bkg[plot_dict['variable']][process].sum().value for process in processes}
    
    sorted_entries = dict(sorted(entries.items(), key=lambda x:x[1], reverse = False))

    # histogram
    bkg_stack = {}
    
    # add histos to stack after rebinning and range setting
    for process in sorted_entries.keys():
        bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][plot_dict['cut'],::style_dict['rebin']]

        # set x range manually
        if style_dict['xlim'] != None:
            xlim = style_dict['xlim']
            xbin_range = np.where((bkg[plot_dict['variable']][process].axes.edges[0] > xlim[0]) & (bkg[plot_dict['variable']][process].axes.edges[0] < xlim[1]))[0]
            bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][ int(xbin_range[0])-1:int(xbin_range[-1]+1) ]

        bkg_stack[process] = bkg[plot_dict['variable']][process]

    hb = hist.Stack.from_dict(bkg_stack)
    
    # Color map for each process
    # cmap = mpl.colormaps['Set3'].colors
    # cmap = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"] # cms-recommended for 6-color scheme
    cmap = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"] # cms-recommended
    
    colors = { 'W+jets': cmap[0],
               'Z+jets': cmap[1],
               'QCD': cmap[2],
               'DY': cmap[3],
               'Top': cmap[4],
               'TTJetsDiLept': cmap[5],
               'Diboson': cmap[6],
               'Triboson': cmap[7],
    }
    
    color_list = [colors[process] for process in sorted_entries.keys()]

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])
    else:
        binwidth = hb[0].axes.widths[0][0]
        
        if style_dict['doDensity']:
            ax.set_ylabel(f'A.U./{binwidth:.3f}')
        else:
            ax.set_ylabel(f'Events/{binwidth:.3f}')

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')
    
    # Plot
    hb.plot(stack=True, yerr=style_dict['doYerr'], density=style_dict['doDensity'], flow=style_dict['flow'], histtype='fill', color=color_list)

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])


def plot_bkg_1d_stacked(ax, bkg_histos, plot_dict, style_dict, processes = 'all'):  

    if processes == 'all':
        #processes = bkg_histos.keys()

        list_cut_index = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=True)
        list_cut_name = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=False)
        
        cut_name = plot_dict['cut']
        
        df = utils.get_bkg_cutflow_df(bkg_histos, 'cutflow_cts').iloc[:-1]
        
        df = df[list_cut_name[list_cut_index.index(cut_name)]]
        
        processes = df.index[df != 0].to_list()
    # if process is given as a list, i.e. ['DY', 'W+jets'], plot only these processes in the list; otherwise, plot all as default
    
    bkg={}
    bkg[plot_dict['variable']] = {process:bkg_histos[process][plot_dict['variable']][{"samp":sum}] for process in processes}
    
    # sort the histograms by the entries and stack
    for process in processes:
        entries = {process: bkg[plot_dict['variable']][process].sum().value for process in processes}
    
    sorted_entries = dict(sorted(entries.items(), key=lambda x:x[1], reverse = False))

    # histogram
    # add histos to stack after rebinning and range setting
    for idx, process in enumerate(sorted_entries.keys()):
        bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][plot_dict['cut'],::style_dict['rebin']]

        # set x range manually
        if style_dict['xlim'] != None:
            xlim = style_dict['xlim']
            xbin_range = np.where((bkg[plot_dict['variable']][process].axes.edges[0] > xlim[0]) & (bkg[plot_dict['variable']][process].axes.edges[0] < xlim[1]))[0]
            bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][ int(xbin_range[0])-1:int(xbin_range[-1]+1) ]
            
        if idx == 0:
            bkg_stack = bkg[plot_dict['variable']][process]
        else:
            bkg_stack += bkg[plot_dict['variable']][process]
    
    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])
    else:
        binwidth = bkg_stack.axes.widths[0][0]
        
        if style_dict['doDensity']:
            ax.set_ylabel(f'A.U./{binwidth:.3f}')
        else:
            ax.set_ylabel(f'Events/{binwidth:.3f}')

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')  
    
    # Plot
    hep.histplot(bkg_stack, yerr=style_dict['doYerr'], density=style_dict['doDensity'], ax=ax, histtype='step', flow=style_dict['flow'], label = style_dict['label'])
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])



def plot_bkg_2D(ax, bkg_histos, plot_dict, style_dict, processes = 'all'):  

    processes_list = processes
    
    if processes == 'all':
        #processes = bkg_histos.keys()

        list_cut_index = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=True)
        list_cut_name = utils.get_bkg_list_of_cuts(bkg_histos, get_cut_idx=False)
        
        cut_name = plot_dict['cut']
        
        df = utils.get_bkg_cutflow_df(bkg_histos, 'cutflow_cts').iloc[:-1]
        
        df = df[list_cut_name[list_cut_index.index(cut_name)]]
        
        processes = df.index[df != 0].to_list()
    # if process is given as a list, i.e. ['DY', 'W+jets'], plot only these processes in the list; otherwise, plot all as default
    
    bkg={}
    bkg[plot_dict['variable']] = {process:bkg_histos[process][plot_dict['variable']][{"samp":sum}] for process in processes}
    
    # sort the histograms by the entries and stack
    for process in processes:
        entries = {process: bkg[plot_dict['variable']][process].sum().value for process in processes}
    
    sorted_entries = dict(sorted(entries.items(), key=lambda x:x[1], reverse = False))

    # histogram
    # add histos to stack after rebinning and range setting
    for idx, process in enumerate(sorted_entries.keys()):
        bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][plot_dict['cut'], ::style_dict['xrebin'], ::style_dict['yrebin']]
    
        if style_dict['xlim'] != None:
            xlim = style_dict['xlim']
            xbin_range = np.where((bkg[plot_dict['variable']][process].axes.edges[0] > xlim[0]) & (bkg[plot_dict['variable']][process].axes.edges[0] < xlim[1]))[0]
            bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][ int(xbin_range[0])-1:int(xbin_range[-1]+1), : ]
        if style_dict['ylim'] != None:
            ylim = style_dict['ylim']
            ybin_range = np.where((bkg[plot_dict['variable']][process].axes.edges[1] > ylim[0]) & (bkg[plot_dict['variable']][process].axes.edges[1] < ylim[1]))[1]
            bkg[plot_dict['variable']][process] = bkg[plot_dict['variable']][process][ :, int(ybin_range[0]):int(ybin_range[-1]+1) ]
            
        if idx == 0:
            bkg_stack = bkg[plot_dict['variable']][process]
        else:
            bkg_stack += bkg[plot_dict['variable']][process]

    # x and y labels
    if style_dict['xlabel'] != None:
        ax.set_xlabel(style_dict['xlabel'])

    if style_dict['ylabel'] != None:
        ax.set_ylabel(style_dict['ylabel'])

    # x,y scale
    if style_dict['doLogx']:
        ax.set_xscale('log')
    if style_dict['doLogy']:
        ax.set_yscale('log')
    
    # Plot
    if style_dict['doLogz']:
        hep.hist2dplot(bkg_stack, flow=style_dict['flow'], norm=mpl.colors.LogNorm(), ax=ax)
    else:
        hep.hist2dplot(bkg_stack, flow=style_dict['flow'], ax=ax)
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
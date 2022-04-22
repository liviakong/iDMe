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

def electron_basics(events,histos,samp):
    eles = events.Electron
    lpt_eles = events.LptElectron
    gen_eles = events.GenPart[np.abs(events.GenPart.ID) == 11]

    ele_kin = histos['ele_kinematics']
    ele_trkHits = histos['ele_trackHits']
    ele_trkQual = histos['ele_trackQual']

    ele_kin.fill(sample=samp,ele_type="Generator",
                        pt=ak.flatten(gen_eles.pt),eta=ak.flatten(gen_eles.eta),phi=ak.flatten(gen_eles.phi))
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
    
    return histos

def ele_genMatching(events,histos,samp):
    eles = events.Electron
    lpt_eles = events.LptElectron
    gen_eles = events.GenPart[np.abs(events.GenPart.ID) == 11]
    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])
    gen_pos = ak.flatten(events.GenPart[events.GenPart.ID == -11])

    match_ele_kin = histos['matched_ele_kinematics']
    match_ele_gen_kin = histos['matched_ele_gen_kinematics']
    match_ele_trkHits = histos['matched_ele_trackHits']
    match_ele_trkQual = histos['matched_ele_trackQual']
    match_ele_gen_disp = histos['matched_ele_gen_displacement']

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
    
    return histos

def genParticles(events,histos,samp):
    gen_eles = events.GenPart[np.abs(events.GenPart.ID) == 11]
    gen_ele_disp = histos['gen_displacement']
    gen_ele_disp.fill(sample=samp,
                        vxy=ak.flatten(gen_eles.vxy),
                        vz=ak.flatten(gen_eles.vz))
    
    return histos

def conversions(events,histos,samp):
    eles = events.Electron
    lpt_eles = events.LptElectron
    convs = events.Conversion
    gen_eles = events.GenPart[np.abs(events.GenPart.ID) == 11]
    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])
    gen_pos = ak.flatten(events.GenPart[events.GenPart.ID == -11])

    match_mtx = histos['match_matrix']
    conv_match_mtx = histos['conv_match_matrix']
    conv_match_mtx_uniq = histos['conv_match_matrix_uniq']

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

    all_eles_dir = ak.concatenate([eles_dir,lpt_eles_dir,conv1_dir,conv2_dir],axis=1)
    all_eles_labels = ak.concatenate([ak.ones_like(eles.eta),2*ak.ones_like(lpt_eles.eta),
                                        3*ak.ones_like(convs.trk1_innerEta),3*ak.ones_like(convs.trk2_innerEta)],axis=1)

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
    ## Matching to most nearby gen object (possible for both conv eles matched to same object)
    d1e_min = ak.fill_none(ak.flatten(d1e[ak.argmin(d1e,axis=1,keepdims=True)]),999)
    d1p_min = ak.fill_none(ak.flatten(d1p[ak.argmin(d1p,axis=1,keepdims=True)]),999)
    d2e_min = ak.fill_none(ak.flatten(d2e[ak.argmin(d2e,axis=1,keepdims=True)]),999)
    d2p_min = ak.fill_none(ak.flatten(d2p[ak.argmin(d2p,axis=1,keepdims=True)]),999)
    matchType_conv1 = ak.values_astype(d1e_min<0.1,int) + 2*ak.values_astype(d1p_min<0.1,int) #0 if no match, 1 if match to e, 2 if match to p, 3 if match to both
    matchType_conv2 = ak.values_astype(d2e_min<0.1,int) + 2*ak.values_astype(d2p_min<0.1,int)
    fullUniqMatch = ((matchType_conv1 == 1) & (matchType_conv2 == 2)) | ((matchType_conv1 == 2) & (matchType_conv2 == 1))
    fullMatch = (matchType_conv1 > 0) & (matchType_conv2 > 0)
    anyMatch = (matchType_conv1 > 0) | (matchType_conv2 > 0)
    conv_match_mtx.fill(sample=samp,
                        match="all",
                        c1_matchType=matchType_conv1,
                        c2_matchType=matchType_conv2)
    conv_match_mtx.fill(sample=samp,
                        match="fullUnique",
                        c1_matchType=matchType_conv1[fullUniqMatch],
                        c2_matchType=matchType_conv2[fullUniqMatch])
    conv_match_mtx.fill(sample=samp,
                        match="full",
                        c1_matchType=matchType_conv1[fullMatch],
                        c2_matchType=matchType_conv2[fullMatch])
    conv_match_mtx.fill(sample=samp,
                        match="any",
                        c1_matchType=matchType_conv1[anyMatch],
                        c2_matchType=matchType_conv2[anyMatch])
    ## Attempting un-ambiguous matches - i.e. "best" matches
    ## If e.g. both conv electrons are best matched to the same object, but one is *better* matched and the 
    ## other could also be matched to the other gen object, just not as well
    both_c1 = (matchType_conv1 == 3)
    both_c2 = (matchType_conv2 == 3)
    # zeroing out the "both" matches
    uniqMatch_conv1 = matchType_conv1*ak.values_astype(~both_c1,int)
    uniqMatch_conv2 = matchType_conv2*ak.values_astype(~both_c2,int)
    # setting to 1 or 2 based on smaller dR, only where there was a "both" match
    uniqMatch_conv1 = uniqMatch_conv1 + ak.values_astype(both_c1,int)*(ak.values_astype(d1e_min < d1p_min,int) + 2*ak.values_astype(d1e_min > d1p_min,int))
    uniqMatch_conv2  = uniqMatch_conv2 + ak.values_astype(both_c2,int)*(ak.values_astype(d2e_min < d2p_min,int) + 2*ak.values_astype(d2e_min > d2p_min,int))
    # Disambiguate some "double" matches -- e.g. conv1 and conv2 both matched to same gen object
    c1c2_e = (uniqMatch_conv1 == uniqMatch_conv2) & (uniqMatch_conv1 == 1)
    c1c2_p = (uniqMatch_conv1 == uniqMatch_conv2) & (uniqMatch_conv1 == 2)
    # zeroing out the "double" matches
    oldMatch_conv1 = 1*uniqMatch_conv1
    oldMatch_conv2 = 1*uniqMatch_conv2
    uniqMatch_conv1 = uniqMatch_conv1*ak.values_astype(~c1c2_e,int)*ak.values_astype(~c1c2_p,int)
    uniqMatch_conv2 = uniqMatch_conv2*ak.values_astype(~c1c2_e,int)*ak.values_astype(~c1c2_p,int)
    # When c1 and c2 are matched to the same object, see if there's a distinct matching possible instead
    switch1_e2p = (c1c2_e) & (d1p_min < 0.1) & ( (d1e_min >= d2e_min) | (d2p_min > 0.1) )
    switch1_p2e = (c1c2_p) & (d1e_min < 0.1) & ( (d1p_min >= d2p_min) | (d2e_min > 0.1) )
    switch2_e2p = (c1c2_e) & (d2p_min < 0.1) & ( (d2e_min > d1e_min) | (d1p_min > 0.1) )
    switch2_p2e = (c1c2_p) & (d2e_min < 0.1) & ( (d2p_min > d1p_min) | (d1e_min > 0.1) )
    no_switch1 = (c1c2_e & (~switch1_e2p)) | (c1c2_p & (~switch1_p2e))
    no_switch2 = (c1c2_e & (~switch2_e2p)) | (c1c2_p & (~switch2_p2e))
    uniqMatch_conv1 = uniqMatch_conv1 + ak.values_astype(switch1_p2e,int) + 2*ak.values_astype(switch1_e2p,int) + oldMatch_conv1*ak.values_astype(no_switch1,int)
    uniqMatch_conv2 = uniqMatch_conv2 + ak.values_astype(switch2_p2e,int) + 2*ak.values_astype(switch2_e2p,int) + oldMatch_conv2*ak.values_astype(no_switch2,int)

    conv_match_mtx_uniq.fill(sample=samp,
                                match="all",
                                c1_matchType=uniqMatch_conv1,
                                c2_matchType=uniqMatch_conv2)
    conv_match_mtx_uniq.fill(sample=samp,
                                match="full",
                                c1_matchType=uniqMatch_conv1[fullMatch],
                                c2_matchType=uniqMatch_conv2[fullMatch])
    conv_match_mtx_uniq.fill(sample=samp,
                                match="any",
                                c1_matchType=uniqMatch_conv1[anyMatch],
                                c2_matchType=uniqMatch_conv2[anyMatch])

    return histos

def recoEE_vertices(events,histos,samp):
    regreg_vtx_vxy = histos["regreg_reco_vertex_vxy"]
    regreg_vtx_stats = histos["regreg_reco_vertex_stats"]
    lowlow_vtx_vxy = histos["lowlow_reco_vertex_vxy"]
    lowlow_vtx_stats = histos["lowlow_reco_vertex_stats"]
    lowreg_vtx_vxy = histos["lowreg_reco_vertex_vxy"]
    lowreg_vtx_stats = histos["lowreg_reco_vertex_stats"]

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
    
    return histos

def recoEE_vertices_genMatch(events,histos,samp):
    reco_vtx_genmatch = histos["reco_vertex_genMatch"]
    nearest_vtx_genmatch = histos["nearest_vertex_genMatch"]
    nearest_vtx_matchType = histos["nearest_vtx_matchType"]

    vtx_regreg = events.EleVertexRegReg
    vtx_lowlow = events.EleVertexLowLow
    vtx_lowreg = events.EleVertexLowReg
    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])

    vtx_regreg_r = ak.zip(
        {
            "x": vtx_regreg.x,
            "y": vtx_regreg.y,
            "z": vtx_regreg.z
        },
        with_name="ThreeVector"
    )

    vtx_lowreg_r = ak.zip(
        {
            "x": vtx_lowreg.x,
            "y": vtx_lowreg.y,
            "z": vtx_lowreg.z
        },
        with_name="ThreeVector"
    )

    vtx_lowlow_r = ak.zip(
        {
            "x": vtx_lowlow.x,
            "y": vtx_lowlow.y,
            "z": vtx_lowlow.z
        },
        with_name="ThreeVector"
    )

    gen_vtx_r = ak.zip(
        {
            "x": gen_ele.x,
            "y": gen_ele.y,
            "z": gen_ele.z
        },
        with_name="ThreeVector"
    )
    
    rr_nearest = ak.argmin((vtx_regreg_r - gen_vtx_r).rho,axis=1,keepdims=True)
    lr_nearest = ak.argmin((vtx_lowreg_r - gen_vtx_r).rho,axis=1,keepdims=True)
    ll_nearest = ak.argmin((vtx_lowlow_r - gen_vtx_r).rho,axis=1,keepdims=True)

    vtx_regreg_nearest = ak.flatten(vtx_regreg[rr_nearest])
    vtx_lowreg_nearest = ak.flatten(vtx_lowreg[lr_nearest])
    vtx_lowlow_nearest = ak.flatten(vtx_lowlow[ll_nearest])
    all_vtx_nearest = ak.concatenate([vtx_regreg[rr_nearest],vtx_lowreg[lr_nearest],vtx_lowlow[ll_nearest]],axis=1)

    dist_rr_nearest = (vtx_regreg_r[rr_nearest] - gen_vtx_r).rho
    dist_lr_nearest = (vtx_lowreg_r[lr_nearest] - gen_vtx_r).rho
    dist_ll_nearest = (vtx_lowlow_r[ll_nearest] - gen_vtx_r).rho

    all_dist_nearest = ak.concatenate([dist_rr_nearest,dist_lr_nearest,dist_ll_nearest],axis=1)
    nearest_ind = ak.argmin(all_dist_nearest,axis=1,keepdims=True)
    global_nearest = ak.flatten(nearest_ind) + 1 # 1 = RR, 2 = LR, 3 = LL
    global_nearest = ak.fill_none(global_nearest,0) # 0 = no vertices in the event!
    has_vtx = global_nearest > 0
    nearest_vtx = ak.flatten(all_vtx_nearest[nearest_ind])
    dist_nearest = ak.flatten(all_dist_nearest[nearest_ind])

    dist_rr_nearest = ak.flatten(dist_rr_nearest)
    dist_lr_nearest = ak.flatten(dist_lr_nearest)
    dist_ll_nearest = ak.flatten(dist_ll_nearest)

    has_rr = ak.count(vtx_regreg.x,axis=1) > 0
    has_lr = ak.count(vtx_lowreg.x,axis=1) > 0
    has_ll = ak.count(vtx_lowlow.x,axis=1) > 0

    # Filling histograms with nearest vertex info split up by type
    reco_vtx_genmatch.fill(sample=samp,
                            type="Reg-Reg",
                            dist=dist_rr_nearest[has_rr],
                            chi2=vtx_regreg_nearest.reduced_chi2[has_rr],
                            dR=vtx_regreg_nearest.dR[has_rr])

    reco_vtx_genmatch.fill(sample=samp,
                            type="Low-Reg",
                            dist=dist_lr_nearest[has_lr],
                            chi2=vtx_lowreg_nearest.reduced_chi2[has_lr],
                            dR=vtx_lowreg_nearest.dR[has_lr])
    
    reco_vtx_genmatch.fill(sample=samp,
                            type="Low-Low",
                            dist=dist_ll_nearest[has_ll],
                            chi2=vtx_lowlow_nearest.reduced_chi2[has_ll],
                            dR=vtx_lowlow_nearest.dR[has_ll])

    # Filling histos with globally nearest vertex
    nearest_vtx_genmatch.fill(sample=samp,
                            type=global_nearest[has_vtx],
                            dist=dist_nearest[has_vtx],
                            chi2=nearest_vtx.reduced_chi2[has_vtx],
                            dR=nearest_vtx.dR[has_vtx])

    nearest_vtx_matchType.fill(sample=samp,
                                type=global_nearest)    

    return histos


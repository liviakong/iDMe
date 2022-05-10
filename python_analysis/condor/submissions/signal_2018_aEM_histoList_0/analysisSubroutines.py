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

def findNearestReco(gen,reco,types,thresh=0.1):
    # construct eta-phi vectors
    gen_r = ak.zip({"x":gen.eta,"y":gen.phi},with_name="TwoVector")
    reco_r = ak.zip({"x":reco.eta,"y":reco.phi},with_name="TwoVector")
    # compute dRs
    dr = (reco_r - gen_r).r  
    # sort by ascending dR
    isort = ak.argsort(dr,axis=1)
    dr = dr[isort]
    types = types[isort]
    # keep only matches with dR < threshold
    match = dr < thresh
    dr = dr[match]
    types = types[match]
    imatch = isort[match] 

    return imatch, types, dr

def findMatches(gen1,gen2,reco,recoType):
    i1, t1, dr1 = findNearestReco(gen1,reco,recoType)
    i2, t2, dr2 = findNearestReco(gen2,reco,recoType)
    return i1, t1, dr1, i2, t2, dr2
    
def findUniqueMatches(gen1,gen2,reco,recoType,allow11=False):
    i1, t1, dr1 = findNearestReco(gen1,reco,recoType)
    i2, t2, dr2 = findNearestReco(gen2,reco,recoType)
    n1 = ak.count(i1,axis=1)
    n2 = ak.count(i2,axis=1)
    
    idx = ak.Array(np.arange(len(reco))) # will help us reconstruct after the fact
    indices1 = [] 
    indices2 = []
    types1 = []
    types2 = []
    drs1 = []
    drs2 = []
    partitions = []

    for i in range(1):
        cutA = (n1 > 0) & (n2 > 0) # remove events where one gen object not matched
        partitions.append(idx[~cutA])
        m1_pad = -1*ak.ones_like(n1[~cutA])[:,np.newaxis]
        i1_pad = ak.concatenate([i1[~cutA],m1_pad],axis=1) # fill -1 where there's no match at all
        i2_pad = ak.concatenate([i2[~cutA],m1_pad],axis=1)
        t1_pad = ak.concatenate([t1[~cutA],m1_pad],axis=1)
        t2_pad = ak.concatenate([t2[~cutA],m1_pad],axis=1)
        dr1_pad = ak.concatenate([dr1[~cutA],m1_pad],axis=1)
        dr2_pad = ak.concatenate([dr2[~cutA],m1_pad],axis=1)
        indices1.append(i1_pad[:,0])
        indices2.append(i2_pad[:,0])
        types1.append(t1_pad[:,0])
        types2.append(t2_pad[:,0])
        drs1.append(dr1_pad[:,0])
        drs2.append(dr2_pad[:,0])
        i1, t1, dr1 = i1[cutA], t1[cutA], dr1[cutA]
        i2, t2, dr2 = i2[cutA], t2[cutA], dr2[cutA]
        idx = idx[cutA]
        n1 = n1[cutA]
        n2 = n2[cutA]
        
        if len(n1) == 0: break
        
        cutB = i1[:,0] == i2[:,0] # remove events where gen objects matched distinctly
        partitions.append(idx[~cutB])
        indices1.append(i1[~cutB][:,0])
        indices2.append(i2[~cutB][:,0])
        types1.append(t1[~cutB][:,0])
        types2.append(t2[~cutB][:,0])
        drs1.append(dr1[~cutB][:,0])
        drs2.append(dr2[~cutB][:,0])
        i1, t1, dr1 = i1[cutB], t1[cutB], dr1[cutB]
        i2, t2, dr2 = i2[cutB], t2[cutB], dr2[cutB]
        idx = idx[cutB]
        n1 = n1[cutB]
        n2 = n2[cutB]

        if len(n1) == 0: break

        cutC = (n1 == 1) & (n2 == 1) # if both matched only to the same thing, then take the match with smaller dR
        cutD = (n1 == 2) & (n2 == 1) 
        cutE = (n1 == 1) & (n2 == 2)
        cutF = (n1 >= 2) & (n2 >= 2)

        if ak.count_nonzero(cutC) > 0: # if both matched only to the same thing, then take the match with smaller dR
            i1c, t1c, dr1c = i1[cutC], t1[cutC], dr1[cutC]
            i2c, t2c, dr2c = i2[cutC], t2[cutC], dr2[cutC]
            if allow11: # option to keep the cases where both have exactly 1 match and it's the same object
                i1c = i1c[:,0]
                t1c = t1c[:,0]
                dr1c = dr1c[:,0]
                i2c = i2c[:,0]
                t2c = t2c[:,0]
                dr2c = dr2c[:,0]
            else:
                d1_nearest = dr1c[:,0]
                d2_nearest = dr2c[:,0]
                i1c = np.where(d1_nearest < d2_nearest,i1c[:,0],-1*ak.ones_like(i1c[:,0]))
                t1c = np.where(d1_nearest < d2_nearest,t1c[:,0],-1*ak.ones_like(t1c[:,0]))
                dr1c = np.where(d1_nearest < d2_nearest,dr1c[:,0],-1*ak.ones_like(dr1c[:,0]))
                i2c = np.where(d2_nearest < d1_nearest,i2c[:,0],-1*ak.ones_like(i2c[:,0]))
                t2c = np.where(d2_nearest < d1_nearest,t2c[:,0],-1*ak.ones_like(t2c[:,0]))
                dr2c = np.where(d2_nearest < d1_nearest,dr2c[:,0],-1*ak.ones_like(dr2c[:,0]))
            partitions.append(idx[cutC])
            indices1.append(i1c)
            indices2.append(i2c)
            types1.append(t1c)
            types2.append(t2c)
            drs1.append(dr1c)
            drs2.append(dr2c)

        if ak.count_nonzero(cutD) > 0: # if there's two options for gen1, take the other option to resolve ambiguity
            i1d, t1d, dr1d = i1[cutD][:,1], t1[cutD][:,1], dr1[cutD][:,1]
            i2d, t2d, dr2d = i2[cutD][:,0], t2[cutD][:,0], dr2[cutD][:,0]
            partitions.append(idx[cutD])
            indices1.append(i1d)
            indices2.append(i2d)
            types1.append(t1d)
            types2.append(t2d)
            drs1.append(dr1d)
            drs2.append(dr2d)

        if ak.count_nonzero(cutE) > 0: # if there's two options for gen2, take the other option to resolve ambiguity
            i2e, t2e, dr2e = i2[cutE][:,1], t2[cutE][:,1], dr2[cutE][:,1]
            i1e, t1e, dr1e = i1[cutE][:,0], t1[cutE][:,0], dr1[cutE][:,0]
            partitions.append(idx[cutE])
            indices1.append(i1e)
            indices2.append(i2e)
            types1.append(t1e)
            types2.append(t2e)
            drs1.append(dr1e)
            drs2.append(dr2e)

        if ak.count_nonzero(cutF) > 0: # if both have multiple options, make the best assignments
            i1f, t1f, dr1f = i1[cutF], t1[cutF], dr1[cutF]
            i2f, t2f, dr2f = i2[cutF], t2[cutF], dr2[cutF]
            d1_nearest = dr1f[:,0] # using slice to preserve structure of the array
            d2_nearest = dr2f[:,0]
            i1f = ak.where(d1_nearest < d2_nearest, i1f[:,0], i1f[:,1])
            t1f = ak.where(d1_nearest < d2_nearest, t1f[:,0], t1f[:,1])
            dr1f = ak.where(d1_nearest < d2_nearest, dr1f[:,0], dr1f[:,1])
            i2f = ak.where(d2_nearest < d1_nearest, i2f[:,0], i2f[:,1])
            t2f = ak.where(d2_nearest < d1_nearest, t2f[:,0], t2f[:,1])
            dr2f = ak.where(d2_nearest < d1_nearest, dr2f[:,0], dr2f[:,1])
            partitions.append(idx[cutF])
            indices1.append(i1f)
            indices2.append(i2f)
            types1.append(t1f)
            types2.append(t2f)
            drs1.append(dr1f)
            drs2.append(dr2f)
    
    all_partitions = ak.concatenate(partitions)
    reorder = ak.argsort(all_partitions)
    type1_out = ak.concatenate(types1)[reorder]
    type2_out = ak.concatenate(types2)[reorder]
    index1_out = ak.concatenate(indices1)[reorder]
    index2_out = ak.concatenate(indices2)[reorder]
    dr1_out = ak.concatenate(drs1)[reorder]
    dr2_out = ak.concatenate(drs2)[reorder]

    return index1_out, type1_out, dr1_out, index2_out, type2_out, dr2_out

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
    match_ele_kin = histos['matched_ele_kinematics']
    match_ele_gen_kin = histos['matched_ele_gen_kinematics']
    match_ele_trkHits = histos['matched_ele_trackHits']
    match_ele_trkQual = histos['matched_ele_trackQual']
    match_ele_gen_disp = histos['matched_ele_gen_displacement']
    match_mtx = histos["ele_match_matrix"]
    nMatch_mtx = histos["ele_nMatch_matrix"]
    match_class = histos["ele_match_class"]

    # Get electron collections
    ele = events.Electron
    lpt_ele = events.LptElectron
    cand_ele = events.EleCand
    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])
    gen_pos = ak.flatten(events.GenPart[events.GenPart.ID == -11])

    # define labels for electron types & concatenate into larger collection
    ele_RL = ak.concatenate([ele,lpt_ele],axis=1)
    type_RL = ak.concatenate([ak.ones_like(ele.pt),2*ak.ones_like(lpt_ele.pt)],axis=1)
    ele_RLC = ak.concatenate([ele_RL,cand_ele],axis=1)
    type_RLC = ak.concatenate([type_RL,3*ak.ones_like(cand_ele.pt)],axis=1)
    
    # matching gen e/p to reco using regular & low-pT electrons
    indE_RL, typeE_RL, drE_RL, indP_RL, typeP_RL, drP_RL = findUniqueMatches(gen_ele,gen_pos,ele_RL,type_RL) # unique matches
    indE_RL_all, typeE_RL_all, drE_RL_all, indP_RL_all, typeP_RL_all, drP_RL_all = findMatches(gen_ele,gen_pos,ele_RL,type_RL) # all matches 
    indE_RL_11, typeE_RL_11, drE_RL_11, indP_RL_11, typeP_RL_11, drP_RL_11 = findUniqueMatches(gen_ele,gen_pos,ele_RL,type_RL,allow11=True) # unique except degenerate where e & p only have 1 match 
    
    # matching gen e/p to reco using regular, low-pT, and candidate electrons
    indE_RLC, typeE_RLC, drE_RLC, indP_RLC, typeP_RLC, drP_RLC = findUniqueMatches(gen_ele,gen_pos,ele_RLC,type_RLC) # unique matches
    indE_RLC_all, typeE_RLC_all, drE_RLC_all, indP_RLC_all, typeP_RLC_all, drP_RLC_all = findMatches(gen_ele,gen_pos,ele_RLC,type_RLC) # all matches 
    indE_RLC_11, typeE_RLC_11, drE_RLC_11, indP_RLC_11, typeP_RLC_11, drP_RLC_11 = findUniqueMatches(gen_ele,gen_pos,ele_RLC,type_RLC,allow11=True) # unique except degenerate where e & p only have 1 match

    # replacing -1s with 0s for the match types
    typeE_RL = np.where(typeE_RL == -1,ak.zeros_like(typeE_RL),typeE_RL)
    typeP_RL = np.where(typeP_RL == -1,ak.zeros_like(typeP_RL),typeP_RL)
    typeE_RL_11 = np.where(typeE_RL_11 == -1,ak.zeros_like(typeE_RL_11),typeE_RL_11)
    typeP_RL_11 = np.where(typeP_RL_11 == -1,ak.zeros_like(typeP_RL_11),typeP_RL_11)
    typeE_RLC = np.where(typeE_RLC == -1,ak.zeros_like(typeE_RLC),typeE_RLC)
    typeP_RLC = np.where(typeP_RLC == -1,ak.zeros_like(typeP_RLC),typeP_RLC)
    typeE_RLC_11 = np.where(typeE_RLC_11 == -1,ak.zeros_like(typeE_RLC_11),typeE_RLC_11)
    typeP_RLC_11 = np.where(typeP_RLC_11 == -1,ak.zeros_like(typeP_RLC_11),typeP_RLC_11)

    # calculating other useful quantities
    nEmatch_RL = ak.count(typeE_RL_all,axis=1)
    nPmatch_RL = ak.count(typeP_RL_all,axis=1)
    nEmatch_RLC = ak.count(typeE_RLC_all,axis=1)
    nPmatch_RLC = ak.count(typeP_RLC_all,axis=1)

    # determining e+e- "match classes" for each event
    zeros = ak.zeros_like(typeE_RL) # generic vector of zeros for each event
    ones = ak.zeros_like(typeE_RL) # generic vector fo ones for each event
    
    noMatch_RL = (typeE_RL_11 == 0) & (typeP_RL_11 == 0)
    oneMatch_RL = ((typeE_RL_11 != 0) & (typeP_RL_11 == 0)) | ((typeE_RL_11 == 0) & (typeP_RL_11 != 0))
    degenMatch_RL = (typeE_RL_11 != 0) & (typeP_RL_11 != 0) & (typeE_RL == typeP_RL)
    uniqMatch_RL = (typeE_RL_11 != 0) & (typeP_RL_11 != 0) & (typeE_RL != typeP_RL)
    matchClass_RL = ak.where(noMatch_RL,zeros,-1*ones)
    matchClass_RL = ak.where(oneMatch_RL,ones,matchClass_RL)
    matchClass_RL = ak.where(degenMatch_RL,2*ones,matchClass_RL)
    matchClass_RL = ak.where(uniqMatch_RL,3*ones,matchClass_RL)

    noMatch_RLC = (typeE_RLC_11 == 0) & (typeP_RLC_11 == 0)
    oneMatch_RLC = ((typeE_RLC_11 != 0) & (typeP_RLC_11 == 0)) | ((typeE_RLC_11 == 0) & (typeP_RLC_11 != 0))
    degenMatch_RLC = (typeE_RLC_11 != 0) & (typeP_RLC_11 != 0) & (typeE_RLC == typeP_RLC)
    uniqMatch_RLC = (typeE_RLC_11 != 0) & (typeP_RLC_11 != 0) & (typeE_RLC != typeP_RLC)
    matchClass_RLC = ak.where(noMatch_RLC,zeros,-1*ones)
    matchClass_RLC = ak.where(oneMatch_RLC,ones,matchClass_RLC)
    matchClass_RLC = ak.where(degenMatch_RLC,2*ones,matchClass_RLC)
    matchClass_RLC = ak.where(uniqMatch_RLC,3*ones,matchClass_RLC)

    # Fill match matrices, reg & low-pT matching
    match_mtx.fill(sample=samp,
                    set="RL",
                    scheme="unique",
                    Etype=typeE_RL,
                    Ptype=typeP_RL)
    match_mtx.fill(sample=samp,
                    set="RL",
                    scheme="nearUnique",
                    Etype=typeE_RL_11,
                    Ptype=typeP_RL_11)
    # Fill match matrices, reg & low-pT & cand matching
    match_mtx.fill(sample=samp,
                    set="RLC",
                    scheme="unique",
                    Etype=typeE_RLC,
                    Ptype=typeP_RLC)
    match_mtx.fill(sample=samp,
                    set="RLC",
                    scheme="nearUnique",
                    Etype=typeE_RLC_11,
                    Ptype=typeP_RLC_11)
    
    # Filling number-of-match matrices
    nMatch_mtx.fill(sample=samp,
                    set="RL",
                    nEmatch=nEmatch_RL,
                    nPmatch=nPmatch_RL)
    nMatch_mtx.fill(sample=samp,
                    set="RLC",
                    nEmatch=nEmatch_RLC,
                    nPmatch=nPmatch_RLC)

    # Filling match class matrices
    match_class.fill(sample=samp,
                    set="RL",
                    matchClass=matchClass_RL)
    match_class.fill(sample=samp,
                    set="RLC",
                    matchClass=matchClass_RLC)

    # Preparing variables and filling matched kinematic histograms
    matched_ele_RL = ak.flatten(ele_RL[indE_RL[indE_RL != -1]])
    matched_pos_RL = ak.flatten(ele_RL[indP_RL[indP_RL != -1]])
    matched_ele_RLC = ak.flatten(ele_RLC[indE_RLC[indE_RLC != -1]])
    matched_pos_RLC = ak.flatten(ele_RLC[indP_RLC[indP_RLC != -1]])
    matches_RL = ak.concatenate([matched_ele_RL,matched_pos_RL])
    matches_RLC = ak.concatenate([matched_ele_RLC,matched_pos_RLC])

    matched_gen_ele_RL = gen_ele[typeE_RL > 0]
    matched_gen_pos_RL = gen_pos[typeP_RL > 0]
    matched_gen_ele_RLC = gen_ele[typeE_RLC > 0]
    matched_gen_pos_RLC = gen_pos[typeP_RLC > 0]
    matches_gen_RL = ak.concatenate([matched_gen_ele_RL,matched_gen_pos_RL])
    matches_gen_RLC = ak.concatenate([matched_gen_ele_RLC,matched_gen_pos_RLC])

    # Filling with reg/lowpT matches
    match_ele_kin.fill(sample=samp,
                        set="RL",
                        pt=matches_RL.pt,
                        eta=matches_RL.eta,
                        phi=matches_RL.phi)
    match_ele_gen_kin.fill(sample=samp,
                            set="RL",
                            pt=matches_gen_RL.pt,
                            eta=matches_gen_RL.eta,
                            phi=matches_gen_RL.phi)
    match_ele_trkHits.fill(sample=samp,
                            set="RL",
                            numTrkHits=matches_RL.numTrackerHits,
                            numPixHits=matches_RL.numPixHits,
                            numStripHits=matches_RL.numStripHits)
    match_ele_trkQual.fill(sample=samp,
                            set="RL",
                            chi2=matches_RL.trkChi2,
                            trkIso=matches_RL.trkIso)
    match_ele_gen_disp.fill(sample=samp,
                            set="RL",
                            vxy=matches_gen_RL.vxy,
                            vz=matches_gen_RL.vz)
    
    # Filling with reg/lowpT/cand matches
    match_ele_kin.fill(sample=samp,
                        set="RLC",
                        pt=matches_RLC.pt,
                        eta=matches_RLC.eta,
                        phi=matches_RLC.phi)
    match_ele_gen_kin.fill(sample=samp,
                            set="RLC",
                            pt=matches_gen_RLC.pt,
                            eta=matches_gen_RLC.eta,
                            phi=matches_gen_RLC.phi)
    match_ele_trkHits.fill(sample=samp,
                            set="RLC",
                            numTrkHits=matches_RLC.numTrackerHits,
                            numPixHits=matches_RLC.numPixHits,
                            numStripHits=matches_RLC.numStripHits)
    match_ele_trkQual.fill(sample=samp,
                            set="RLC",
                            chi2=matches_RLC.trkChi2,
                            trkIso=matches_RLC.trkIso)
    match_ele_gen_disp.fill(sample=samp,
                            set="RLC",
                            vxy=matches_gen_RLC.vxy,
                            vz=matches_gen_RLC.vz)
    
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
    rr_vtx_vxy = histos["RRvtx_vxy"]
    rr_vtx_stats = histos["RRvtx_stats"]
    ll_vtx_vxy = histos["LLvtx_vxy"]
    ll_vtx_stats = histos["LLvtx_stats"]
    lr_vtx_vxy = histos["LRvtx_vxy"]
    lr_vtx_stats = histos["LRvtx_stats"]

    ################ Reco Electron Vertex Histograms ################
    vtx_rr = events.RRvtx
    vtx_ll = events.LLvtx
    vtx_lr = events.LRvtx
    # Selecting only the entries with a valid vertex, i.e. vxy > 0
    vtx_rr = vtx_rr[vtx_rr.vxy > 0]
    vtx_ll = vtx_ll[vtx_ll.vxy > 0]
    vtx_lr = vtx_lr[vtx_lr.vxy > 0]
    # Filters for selecting opposite charge dielectron vertices
    opp_rr = vtx_rr.sign == -1
    opp_ll = vtx_ll.sign == -1
    opp_lr = vtx_lr.sign == -1
    # reg-reg, i.e. two regular electrons matched to a vertex
    rr_vtx_vxy.fill(sample=samp,
                        sign="all",
                        vxy=ak.flatten(vtx_rr.vxy),
                        sigma_vxy=ak.flatten(vtx_rr.sigmavxy))
    rr_vtx_vxy.fill(sample=samp,
                        sign="opp",
                        vxy=ak.flatten(vtx_rr[opp_rr].vxy),
                        sigma_vxy=ak.flatten(vtx_rr[opp_rr].sigmavxy))
    rr_vtx_stats.fill(sample=samp,
                            sign="all",
                            chi2=ak.flatten(vtx_rr.reduced_chi2),
                            signif=ak.flatten(vtx_rr.vxy/vtx_rr.sigmavxy))
    rr_vtx_stats.fill(sample=samp,
                            sign="opp",
                            chi2=ak.flatten(vtx_rr[opp_rr].reduced_chi2),
                            signif=ak.flatten(vtx_rr[opp_rr].vxy/vtx_rr[opp_rr].sigmavxy))
    # low-low, i.e. two low-pT electrons matched to a vertex
    ll_vtx_vxy.fill(sample=samp,
                        sign="all",
                        vxy=ak.flatten(vtx_ll.vxy),
                        sigma_vxy=ak.flatten(vtx_ll.sigmavxy))
    ll_vtx_vxy.fill(sample=samp,
                        sign="opp",
                        vxy=ak.flatten(vtx_ll[opp_ll].vxy),
                        sigma_vxy=ak.flatten(vtx_ll[opp_ll].sigmavxy))
    ll_vtx_stats.fill(sample=samp,
                            sign="all",
                            chi2=ak.flatten(vtx_ll.reduced_chi2),
                            signif=ak.flatten(vtx_ll.vxy/vtx_ll.sigmavxy))
    ll_vtx_stats.fill(sample=samp,
                            sign="opp",
                            chi2=ak.flatten(vtx_ll[opp_ll].reduced_chi2),
                            signif=ak.flatten(vtx_ll[opp_ll].vxy/vtx_ll[opp_ll].sigmavxy))
    # low-reg, i.e. one low-pT and one default electron matched to a vertex
    lr_vtx_vxy.fill(sample=samp,
                        sign="all",
                        vxy=ak.flatten(vtx_lr.vxy),
                        sigma_vxy=ak.flatten(vtx_lr.sigmavxy))
    lr_vtx_vxy.fill(sample=samp,
                        sign="opp",
                        vxy=ak.flatten(vtx_lr[opp_lr].vxy),
                        sigma_vxy=ak.flatten(vtx_lr[opp_lr].sigmavxy))
    lr_vtx_stats.fill(sample=samp,
                            sign="all",
                            chi2=ak.flatten(vtx_lr.reduced_chi2),
                            signif=ak.flatten(vtx_lr.vxy/vtx_lr.sigmavxy))
    lr_vtx_stats.fill(sample=samp,
                            sign="opp",
                            chi2=ak.flatten(vtx_lr[opp_lr].reduced_chi2),
                            signif=ak.flatten(vtx_lr[opp_lr].vxy/vtx_lr[opp_lr].sigmavxy))
    
    return histos

def recoEE_vertices_genMatchDr(events,histos,samp):
    reco_vtx_genmatch = histos["vtx_genMatchByDr"]
    reco_vtx_kin_genmatch = histos["vtx_kin_genMatchByDr"]
    reco_vtx_disp_genmatch = histos["vtx_disp_genMatchByDr"]
    nearest_vtx_genmatch = histos["vtx_nearestMatchByDr"]
    nearest_vtx_kin_genmatch = histos["vtx_kin_nearestMatchByDr"]
    nearest_vtx_disp_genmatch = histos["vtx_disp_nearestMatchByDr"]
    nearest_vtx_matchType = histos["vtx_nearestMatchTypeByDr"]

    vtx_rr = events.RRvtx
    vtx_ll = events.LLvtx
    vtx_lr = events.LRvtx
    vtx_rc = events.RCvtx
    vtx_lc = events.LCvtx

    n_rr = ak.count(vtx_rr.vx,axis=1)
    n_ll = ak.count(vtx_ll.vx,axis=1)
    n_lr = ak.count(vtx_lr.vx,axis=1)
    n_rc = ak.count(vtx_rc.vx,axis=1)
    n_lc = ak.count(vtx_lc.vx,axis=1)

    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])

    vtx_rr_r = ak.zip(
        {
            "x": vtx_rr.vx,
            "y": vtx_rr.vy,
            "z": vtx_rr.vz
        },
        with_name="ThreeVector"
    )

    vtx_lr_r = ak.zip(
        {
            "x": vtx_lr.vx,
            "y": vtx_lr.vy,
            "z": vtx_lr.vz
        },
        with_name="ThreeVector"
    )

    vtx_ll_r = ak.zip(
        {
            "x": vtx_ll.vx,
            "y": vtx_ll.vy,
            "z": vtx_ll.vz
        },
        with_name="ThreeVector"
    )

    vtx_rc_r = ak.zip(
        {
            "x": vtx_rc.vx,
            "y": vtx_rc.vy,
            "z": vtx_rc.vz
        },
        with_name="ThreeVector"
    )

    vtx_lc_r = ak.zip(
        {
            "x": vtx_lc.vx,
            "y": vtx_lc.vy,
            "z": vtx_lc.vz
        },
        with_name="ThreeVector"
    )

    gen_vtx_r = ak.zip(
        {
            "x": gen_ele.vx,
            "y": gen_ele.vy,
            "z": gen_ele.vz
        },
        with_name="ThreeVector"
    )
    
    rr_nearest = ak.argmin((vtx_rr_r - gen_vtx_r).rho,axis=1,keepdims=True)
    lr_nearest = ak.argmin((vtx_lr_r - gen_vtx_r).rho,axis=1,keepdims=True)
    ll_nearest = ak.argmin((vtx_ll_r - gen_vtx_r).rho,axis=1,keepdims=True)

    vtx_rr_nearest = ak.flatten(vtx_rr[rr_nearest])
    vtx_lr_nearest = ak.flatten(vtx_lr[lr_nearest])
    vtx_ll_nearest = ak.flatten(vtx_ll[ll_nearest])
    all_vtx_nearest = ak.concatenate([vtx_rr[rr_nearest],vtx_lr[lr_nearest],vtx_ll[ll_nearest]],axis=1)

    dist_rr_nearest = (vtx_rr_r[rr_nearest] - gen_vtx_r).rho
    dist_lr_nearest = (vtx_lr_r[lr_nearest] - gen_vtx_r).rho
    dist_ll_nearest = (vtx_ll_r[ll_nearest] - gen_vtx_r).rho

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

    has_rr = ak.count(vtx_rr.vx,axis=1) > 0
    has_lr = ak.count(vtx_lr.vx,axis=1) > 0
    has_ll = ak.count(vtx_ll.vx,axis=1) > 0

    # Filling histograms with nearest vertex info split up by type
    reco_vtx_genmatch.fill(sample=samp,
                            type="Reg-Reg",
                            dist=dist_rr_nearest[has_rr],
                            chi2=vtx_rr_nearest.reduced_chi2[has_rr],
                            dR=vtx_rr_nearest.dR[has_rr])
    reco_vtx_genmatch.fill(sample=samp,
                            type="Low-Reg",
                            dist=dist_lr_nearest[has_lr],
                            chi2=vtx_lr_nearest.reduced_chi2[has_lr],
                            dR=vtx_lr_nearest.dR[has_lr])
    reco_vtx_genmatch.fill(sample=samp,
                            type="Low-Low",
                            dist=dist_ll_nearest[has_ll],
                            chi2=vtx_ll_nearest.reduced_chi2[has_ll],
                            dR=vtx_ll_nearest.dR[has_ll])

    reco_vtx_kin_genmatch.fill(sample=samp,
                                type="Reg-Reg",
                                mass=vtx_rr_nearest.m[has_rr],
                                energy=vtx_rr_nearest.energy[has_rr])
    reco_vtx_kin_genmatch.fill(sample=samp,
                                type="Low-Reg",
                                mass=vtx_lr_nearest.m[has_lr],
                                energy=vtx_lr_nearest.energy[has_lr])
    reco_vtx_kin_genmatch.fill(sample=samp,
                                type="Low-Low",
                                mass=vtx_ll_nearest.m[has_ll],
                                energy=vtx_ll_nearest.energy[has_ll])

    reco_vtx_disp_genmatch.fill(sample=samp,
                                type="Reg-Reg",
                                vxy=vtx_rr_nearest.vxy[has_rr],
                                vxy_signif=(vtx_rr_nearest.vxy/vtx_rr_nearest.sigmavxy)[has_rr])
    reco_vtx_disp_genmatch.fill(sample=samp,
                                type="Low-Reg",
                                vxy=vtx_lr_nearest.vxy[has_lr],
                                vxy_signif=(vtx_lr_nearest.vxy/vtx_lr_nearest.sigmavxy)[has_lr])
    reco_vtx_disp_genmatch.fill(sample=samp,
                                type="Low-Low",
                                vxy=vtx_ll_nearest.vxy[has_ll],
                                vxy_signif=(vtx_ll_nearest.vxy/vtx_ll_nearest.sigmavxy)[has_ll])

    # Filling histos with globally nearest vertex
    nearest_vtx_genmatch.fill(sample=samp,
                            type=global_nearest[has_vtx],
                            dist=dist_nearest[has_vtx],
                            chi2=nearest_vtx.reduced_chi2[has_vtx],
                            dR=nearest_vtx.dR[has_vtx])
    nearest_vtx_kin_genmatch.fill(sample=samp,
                                  type=global_nearest[has_vtx],
                                  mass=nearest_vtx.m[has_vtx],
                                  energy=nearest_vtx.energy[has_vtx])
    nearest_vtx_disp_genmatch.fill(sample=samp,
                                   type=global_nearest[has_vtx],
                                   vxy=nearest_vtx.vxy[has_vtx],
                                   vxy_signif=(nearest_vtx.vxy/nearest_vtx.sigmavxy)[has_vtx])

    nearest_vtx_matchType.fill(sample=samp,
                                type=global_nearest)


    # Looking at displaced dilepton vertices found by standalone algorithms
    h_dEE_kin = histos["dispEE_vtx_kinematics"]
    h_dEE_qual = histos["dispEE_vtx_qual"]
    h_dEE_eff = histos["dispEE_eff"]

    dEE = events.EECand
    dEE_vxy = np.sqrt(dEE.vx**2 + dEE.vy**2)
    ndEE = ak.count(dEE.dR,axis=1)
    has_dEE = ndEE > 0

    eles_r = ak.zip({"x":events.Electron.eta,"y":events.Electron.phi},with_name="TwoVector")
    lpt_eles_r = ak.zip({"x":events.LptElectron.eta,"y":events.LptElectron.phi},with_name="TwoVector")
    gen_ele = events.GenPart[events.GenPart.ID == 11]
    gen_pos = events.GenPart[events.GenPart.ID == -11]
    gen_ele_r = ak.flatten(ak.zip({"x":gen_ele.eta,"y":gen_ele.phi},with_name="TwoVector"))
    gen_pos_r = ak.flatten(ak.zip({"x":gen_pos.eta,"y":gen_pos.phi},with_name="TwoVector"))

    mindr_ee = ak.min((eles_r - gen_ele_r).r,axis=1)
    mindr_ep = ak.min((eles_r - gen_pos_r).r,axis=1)
    mindr_lee = ak.min((lpt_eles_r - gen_ele_r).r,axis=1)
    mindr_lep = ak.min((lpt_eles_r - gen_pos_r).r,axis=1)

    # Masks for denominator categories
    noRegVtx = (n_rr == 0) & (n_lr == 0) & (n_ll == 0)
    noGenMatches = ak.fill_none((mindr_ee > 0.1) & (mindr_ep > 0.1) & (mindr_lee > 0.1) & (mindr_lep > 0.1),True)

    # inclusive
    h_dEE_eff.fill(sample=samp,
                    denom="all",
                    ndEE=ndEE)
    h_dEE_kin.fill(sample=samp,
                    denom="all",
                    mass=ak.flatten(dEE[has_dEE].mass),
                    leadPt=ak.flatten(dEE[has_dEE].leadingPt),
                    dR=ak.flatten(dEE[has_dEE].dR))
    h_dEE_qual.fill(sample=samp,
                    denom="all",
                    vxy=ak.flatten(dEE_vxy[has_dEE]),
                    dxy=ak.flatten(dEE[has_dEE].trackDxy_PV),
                    chi2=ak.flatten(dEE[has_dEE].normalizedChi2))

    # events with no RR, LR, or LL vertices
    h_dEE_eff.fill(sample=samp,
                    denom="noRegVtx",
                    ndEE=ndEE[noRegVtx])
    h_dEE_kin.fill(sample=samp,
                    denom="noRegVtx",
                    mass=ak.flatten(dEE[has_dEE & noRegVtx].mass),
                    leadPt=ak.flatten(dEE[has_dEE & noRegVtx].leadingPt),
                    dR=ak.flatten(dEE[has_dEE & noRegVtx].dR))
    h_dEE_qual.fill(sample=samp,
                    denom="noRegVtx",
                    vxy=ak.flatten(dEE_vxy[has_dEE & noRegVtx]),
                    dxy=ak.flatten(dEE[has_dEE & noRegVtx].trackDxy_PV),
                    chi2=ak.flatten(dEE[has_dEE & noRegVtx].normalizedChi2))

    # Events with no gen-matched electrons (reg or low-pt)
    h_dEE_eff.fill(sample=samp,
                    denom="noGenMatchEles",
                    ndEE=ndEE[noGenMatches])
    h_dEE_kin.fill(sample=samp,
                    denom="noGenMatchEles",
                    mass=ak.flatten(dEE[has_dEE & noGenMatches].mass),
                    leadPt=ak.flatten(dEE[has_dEE & noGenMatches].leadingPt),
                    dR=ak.flatten(dEE[has_dEE & noGenMatches].dR))
    h_dEE_qual.fill(sample=samp,
                    denom="noGenMatchEles",
                    vxy=ak.flatten(dEE_vxy[has_dEE & noGenMatches]),
                    dxy=ak.flatten(dEE[has_dEE & noGenMatches].trackDxy_PV),
                    chi2=ak.flatten(dEE[has_dEE & noGenMatches].normalizedChi2))

    return histos

def recoEE_vertices_genMatchByEles(events,histos,samp):
    # Reco vertices
    vtx_rr = events.RRvtx
    vtx_ll = events.LLvtx
    vtx_lr = events.LRvtx
    vtx_rc = events.RCvtx
    vtx_lc = events.LCvtx

    ele = events.Electron
    lpt_ele = events.LptElectron
    cand_ele = events.EleCand
    gen_ele = ak.flatten(events.GenPart[events.GenPart.ID == 11])
    gen_pos = ak.flatten(events.GenPart[events.GenPart.ID == -11])
    # define labels for electron types & concatenate into larger collection
    ele_RL = ak.concatenate([ele,lpt_ele],axis=1)
    type_RL = ak.concatenate([ak.ones_like(ele.pt),2*ak.ones_like(lpt_ele.pt)],axis=1)
    ele_RLC = ak.concatenate([ele_RL,cand_ele],axis=1)
    type_RLC = ak.concatenate([type_RL,3*ak.ones_like(cand_ele.pt)],axis=1)
    # matching gen e/p to reco using regular & low-pT electrons
    indE_RL, typeE_RL, drE_RL, indP_RL, typeP_RL, drP_RL = findUniqueMatches(gen_ele,gen_pos,ele_RL,type_RL) # unique matches
    # matching gen e/p to reco using regular, low-pT, and candidate electrons
    indE_RLC, typeE_RLC, drE_RLC, indP_RLC, typeP_RLC, drP_RLC = findUniqueMatches(gen_ele,gen_pos,ele_RLC,type_RLC) # unique matches

    # replacing -1s with 0s for the match types
    typeE_RL = np.where(typeE_RL == -1,ak.zeros_like(typeE_RL),typeE_RL)
    typeP_RL = np.where(typeP_RL == -1,ak.zeros_like(typeP_RL),typeP_RL)
    typeE_RLC = np.where(typeE_RLC == -1,ak.zeros_like(typeE_RLC),typeE_RLC)
    typeP_RLC = np.where(typeP_RLC == -1,ak.zeros_like(typeP_RLC),typeP_RLC)

    # masks for different vertex types, low/reg only
    matchRR_RL = (typeE_RL == 1) & (typeP_RL == 1)
    matchLR_RL = ((typeE_RL == 1) & (typeP_RL == 2)) | ((typeE_RL == 2) & (typeP_RL == 1))
    matchLL_RL = (typeE_RL == 2) & (typeP_RL == 2)

    # masks for different vertex types, low/reg/cand
    matchRR_RLC = (typeE_RLC == 1) & (typeP_RLC == 1)
    matchRC_RLC = ((typeE_RLC == 1) & (typeP_RLC == 3)) | ((typeE_RLC == 3) & (typeP_RLC == 1))
    matchLR_RLC = ((typeE_RLC == 1) & (typeP_RLC == 2)) | ((typeE_RLC == 2) & (typeP_RLC == 1))
    matchLC_RLC = ((typeE_RLC == 2) & (typeP_RLC == 3)) | ((typeE_RLC == 3) & (typeP_RLC == 2))
    matchLL_RLC = (typeE_RLC == 2) & (typeP_RLC == 2)
    matchCC_RLC = (typeE_RLC == 3) & (typeP_RLC == 3)

    # Retrieving the matched vertices for each event
    n_reg_eles = ak.count(ele.pt,axis=1) # for retrieving "local" indices within individual collections
    n_lpt_eles = ak.count(lpt_ele.pt,axis=1)
    n_cand_eles = ak.count(cand_ele.pt,axis=1)

    # reg/lowpT only
    RRvtx_RL = vtx_rr[matchRR_RL]
    RRvtx_mask_RL = ((indE_RL[matchRR_RL] == RRvtx_RL.idx1) & (indP_RL[matchRR_RL] == RRvtx_RL.idx2)) |\
                    ((indE_RL[matchRR_RL] == RRvtx_RL.idx2) & (indP_RL[matchRR_RL] == RRvtx_RL.idx1))
    RRvtx_RL = ak.flatten(RRvtx_RL[RRvtx_mask_RL])

    LRvtx_RL = vtx_lr[matchLR_RL]
    LRvtx_mask_RL = (((indE_RL[matchLR_RL] - n_reg_eles[matchLR_RL]) == LRvtx_RL.idx1) & (indP_RL[matchLR_RL] == LRvtx_RL.idx2)) |\
                    ((indE_RL[matchLR_RL] == LRvtx_RL.idx2) & ((indP_RL[matchLR_RL] - n_reg_eles[matchLR_RL]) == LRvtx_RL.idx1))
    LRvtx_RL = ak.flatten(LRvtx_RL[LRvtx_mask_RL])

    LLvtx_RL = vtx_ll[matchLL_RL]
    LLvtx_mask_RL = ((indE_RL[matchLL_RL] == LLvtx_RL.idx1) & (indP_RL[matchLL_RL] == LLvtx_RL.idx2)) |\
                    ((indE_RL[matchLL_RL] == LLvtx_RL.idx2) & (indP_RL[matchLL_RL] == LLvtx_RL.idx1))
    LLvtx_RL = ak.flatten(LLvtx_RL[LLvtx_mask_RL])

    # reg,lowpT,cand
    RRvtx_RLC = vtx_rr[matchRR_RLC]
    RRvtx_mask_RLC = ((indE_RLC[matchRR_RLC] == RRvtx_RLC.idx1) & (indP_RLC[matchRR_RLC] == RRvtx_RLC.idx2)) |\
                     ((indE_RLC[matchRR_RLC] == RRvtx_RLC.idx2) & (indP_RLC[matchRR_RLC] == RRvtx_RLC.idx1))
    RRvtx_RLC = ak.flatten(RRvtx_RLC[RRvtx_mask_RLC])

    RCvtx_RLC = vtx_rc[matchRC_RLC]
    RCvtx_mask_RLC = ((indE_RLC[matchRC_RLC] == RCvtx_RLC.idx1) & ((indP_RLC[matchRC_RLC] - n_reg_eles[matchRC_RLC] - n_lpt_eles[matchRC_RLC]) == RCvtx_RLC.idx2)) |\
                     (((indE_RLC[matchRC_RLC] - n_reg_eles[matchRC_RLC] - n_lpt_eles[matchRC_RLC]) == RCvtx_RLC.idx2) & (indP_RLC[matchRC_RLC] == RCvtx_RLC.idx1))
    RCvtx_RLC = ak.flatten(RCvtx_RLC[RCvtx_mask_RLC])

    LCvtx_RLC = vtx_lc[matchLC_RLC]
    LCvtx_mask_RLC = (((indE_RLC[matchLC_RLC] - n_reg_eles[matchLC_RLC]) == LCvtx_RLC.idx1) & ((indP_RLC[matchLC_RLC] - n_reg_eles[matchLC_RLC] - n_lpt_eles[matchLC_RLC]) == LCvtx_RLC.idx2)) |\
                     (((indE_RLC[matchLC_RLC] - n_reg_eles[matchLC_RLC] - n_lpt_eles[matchLC_RLC]) == LCvtx_RLC.idx2) & ((indP_RLC[matchLC_RLC] - n_reg_eles[matchLC_RLC]) == LCvtx_RLC.idx1))
    LCvtx_RLC = ak.flatten(LCvtx_RLC[LCvtx_mask_RLC])

    LRvtx_RLC = vtx_lr[matchLR_RLC]
    LRvtx_mask_RLC = (((indE_RLC[matchLR_RLC] - n_reg_eles[matchLR_RLC]) == LRvtx_RLC.idx1) & (indP_RLC[matchLR_RLC] == LRvtx_RLC.idx2)) |\
                     ((indE_RLC[matchLR_RLC] == LRvtx_RLC.idx2) & ((indP_RLC[matchLR_RLC] - n_reg_eles[matchLR_RLC]) == LRvtx_RLC.idx1))
    LRvtx_RLC = ak.flatten(LRvtx_RLC[LRvtx_mask_RLC])

    LLvtx_RLC = vtx_ll[matchLL_RLC]
    LLvtx_mask_RLC = ((indE_RLC[matchLL_RLC] == LLvtx_RLC.idx1) & (indP_RLC[matchLL_RLC] == LLvtx_RLC.idx2)) |\
                     ((indE_RLC[matchLL_RLC] == LLvtx_RLC.idx2) & (indP_RLC[matchLL_RLC] == LLvtx_RLC.idx1))
    LLvtx_RLC = ak.flatten(LLvtx_RLC[LLvtx_mask_RLC])

    # Filling histograms, finally!
    h_vtx = histos["vtx_genMatchByEle"]
    h_vtx_kin = histos["vtx_kin_genMatchByEle"]
    h_vtx_disp = histos["vtx_disp_genMatchByEle"]

    # filling with the reco/lowpT set first
    h_vtx.fill(sample=samp,
               type="Reg-Reg",
               set="RL",
               chi2=RRvtx_RL.reduced_chi2,
               dR=RRvtx_RL.dR)
    h_vtx.fill(sample=samp,
               type="Low-Reg",
               set="RL",
               chi2=LRvtx_RL.reduced_chi2,
               dR=LRvtx_RL.dR)
    h_vtx.fill(sample=samp,
               type="Low-Low",
               set="RL",
               chi2=LLvtx_RL.reduced_chi2,
               dR=LLvtx_RL.dR)
    
    h_vtx_kin.fill(sample=samp,
                   type="Reg-Reg",
                   set="RL",
                   mass=RRvtx_RL.m,
                   energy=RRvtx_RL.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Low-Reg",
                   set="RL",
                   mass=LRvtx_RL.m,
                   energy=LRvtx_RL.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Low-Low",
                   set="RL",
                   mass=LLvtx_RL.m,
                   energy=LLvtx_RL.energy)
    
    h_vtx_disp.fill(sample=samp,
                    type="Reg-Reg",
                    set="RL",
                    vxy=RRvtx_RL.vxy,
                    vxy_signif=(RRvtx_RL.vxy/RRvtx_RL.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Low-Reg",
                    set="RL",
                    vxy=LRvtx_RL.vxy,
                    vxy_signif=(LRvtx_RL.vxy/LRvtx_RL.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Low-Low",
                    set="RL",
                    vxy=LLvtx_RL.vxy,
                    vxy_signif=(LLvtx_RL.vxy/LLvtx_RL.sigmavxy))

    # Filling with reg/lowpT/cand set
    h_vtx.fill(sample=samp,
               type="Reg-Reg",
               set="RLC",
               chi2=RRvtx_RLC.reduced_chi2,
               dR=RRvtx_RLC.dR)
    h_vtx.fill(sample=samp,
               type="Reg-Cand",
               set="RLC",
               chi2=RCvtx_RLC.reduced_chi2,
               dR=RCvtx_RLC.dR)
    h_vtx.fill(sample=samp,
               type="Low-Reg",
               set="RLC",
               chi2=LRvtx_RLC.reduced_chi2,
               dR=LRvtx_RLC.dR)
    h_vtx.fill(sample=samp,
               type="Low-Cand",
               set="RLC",
               chi2=LCvtx_RLC.reduced_chi2,
               dR=LCvtx_RLC.dR)
    h_vtx.fill(sample=samp,
               type="Low-Low",
               set="RLC",
               chi2=LLvtx_RLC.reduced_chi2,
               dR=LLvtx_RLC.dR)
    
    h_vtx_kin.fill(sample=samp,
                   type="Reg-Reg",
                   set="RLC",
                   mass=RRvtx_RLC.m,
                   energy=RRvtx_RLC.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Reg-Cand",
                   set="RLC",
                   mass=RCvtx_RLC.m,
                   energy=RCvtx_RLC.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Low-Reg",
                   set="RLC",
                   mass=LRvtx_RLC.m,
                   energy=LRvtx_RLC.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Low-Cand",
                   set="RLC",
                   mass=LCvtx_RLC.m,
                   energy=LCvtx_RLC.energy)
    h_vtx_kin.fill(sample=samp,
                   type="Low-Low",
                   set="RLC",
                   mass=LLvtx_RLC.m,
                   energy=LLvtx_RLC.energy)
    
    h_vtx_disp.fill(sample=samp,
                    type="Reg-Reg",
                    set="RLC",
                    vxy=RRvtx_RLC.vxy,
                    vxy_signif=(RRvtx_RLC.vxy/RRvtx_RLC.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Reg-Cand",
                    set="RLC",
                    vxy=RCvtx_RLC.vxy,
                    vxy_signif=(RCvtx_RLC.vxy/RCvtx_RLC.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Low-Reg",
                    set="RLC",
                    vxy=LRvtx_RLC.vxy,
                    vxy_signif=(LRvtx_RLC.vxy/LRvtx_RLC.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Low-Cand",
                    set="RLC",
                    vxy=LCvtx_RLC.vxy,
                    vxy_signif=(LCvtx_RLC.vxy/LCvtx_RLC.sigmavxy))
    h_vtx_disp.fill(sample=samp,
                    type="Low-Low",
                    set="RLC",
                    vxy=LLvtx_RLC.vxy,
                    vxy_signif=(LLvtx_RLC.vxy/LLvtx_RLC.sigmavxy))
    
    return histos
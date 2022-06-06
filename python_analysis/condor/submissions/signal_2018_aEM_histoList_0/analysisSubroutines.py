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

def findUniqueMatches(ele_matches,best_ele_match,pos_matches,best_pos_match,useCands=False,allow11=False):
    if not useCands:
        ele_matches = ele_matches[ele_matches.typ != 3]
        pos_matches = pos_matches[pos_matches.typ != 3]
    numE = ak.count(ele_matches.typ,axis=1)
    numP = ak.count(pos_matches.typ,axis=1)
    idx = ak.Array(np.arange(len(best_ele_match)))

    partitions = []
    ele = []
    pos = []

    for i in range(1):
        # removing events where one or more objects has no match
        cutA = (numE > 0) & (numP > 0)
        ele.append(best_ele_match[~cutA])
        pos.append(best_pos_match[~cutA])
        partitions.append(idx[~cutA])
        ele_matches = ele_matches[cutA]
        best_ele_match = best_ele_match[cutA]
        pos_matches = pos_matches[cutA]
        best_pos_match = best_pos_match[cutA]
        idx = idx[cutA]
        numE = numE[cutA]
        numP = numP[cutA]

        if len(numE) == 0:
            break

        # remove events with distinct match
        cutB = (best_ele_match.typ == best_pos_match.typ) & (best_ele_match.ind == best_pos_match.ind)
        ele.append(best_ele_match[~cutB])
        pos.append(best_pos_match[~cutB])
        partitions.append(idx[~cutB])
        ele_matches = ele_matches[cutB]
        best_ele_match = best_ele_match[cutB]
        pos_matches = pos_matches[cutB]
        best_pos_match = best_pos_match[cutB]
        idx = idx[cutB]
        numE = numE[cutB]
        numP = numP[cutB]

        if len(numP) == 0:
            break

        # sorting remaining matches by dR
        sort_ele = ak.argsort(ele_matches.dr,axis=1)
        sort_pos = ak.argsort(pos_matches.dr,axis=1)
        ele_matches = ele_matches[sort_ele]
        pos_matches = pos_matches[sort_pos]
        
        cutC = (numE == 1) & (numP == 1)
        cutD = (numE >= 2) & (numP == 1)
        cutE = (numE == 1) & (numP >= 2)
        cutF = (numE >= 2) & (numP >= 2)

        if ak.count_nonzero(cutC) > 0:
            ele_best = best_ele_match[cutC]
            pos_best = best_pos_match[cutC]
            if not allow11:
                ele_best = ak.where(ele_best.dr < pos_best.dr,ele_best,ak.full_like(ele_best,0))
                pos_best = ak.where(pos_best.dr < ele_best.dr,pos_best,ak.full_like(pos_best,0))
            ele.append(ele_best)
            pos.append(pos_best)
            partitions.append(idx[cutC])
        if ak.count_nonzero(cutD) > 0:
            ele_best = ele_matches[cutD][:,1]
            pos_best = best_pos_match[cutD]
            ele.append(ele_best)
            pos.append(pos_best)
            partitions.append(idx[cutD])
        if ak.count_nonzero(cutE) > 0:
            ele_best = best_ele_match[cutE]
            pos_best = pos_matches[cutE][:,1]
            ele.append(ele_best)
            pos.append(pos_best)
            partitions.append(idx[cutE])
        if ak.count_nonzero(cutF) > 0:
            ele_best = ak.where(best_ele_match[cutF].dr < best_pos_match[cutF].dr,ele_matches[cutF][:,0],ele_matches[cutF][:,1])
            pos_best = ak.where(best_pos_match[cutF].dr < best_ele_match[cutF].dr,pos_matches[cutF][:,0],pos_matches[cutF][:,1])
            ele.append(ele_best)
            pos.append(pos_best)
            partitions.append(idx[cutF])
    
    all_partitions = ak.concatenate(partitions)
    reorder = ak.argsort(all_partitions)
    ele_out = ak.concatenate(ele)[reorder]
    pos_out = ak.concatenate(pos)[reorder]

    return ele_out, pos_out

def findUniqueMatches_old(gen1,gen2,reco,recoType,allow11=False):
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
        cutD = (n1 >= 2) & (n2 == 1) 
        cutE = (n1 == 1) & (n2 >= 2)
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
    cand_eles = events.EleCand

    ele_kin = histos['ele_kinematics']
    ele_trkHits = histos['ele_trackHits']
    ele_trkQual = histos['ele_trackQual']

    ele_kin.fill(sample=samp,
                ele_type="Default",
                pt=ak.flatten(eles.pt),
                eta=ak.flatten(eles.eta),
                phi=ak.flatten(eles.phi))
    ele_kin.fill(sample=samp,
                ele_type="Low pT",
                pt=ak.flatten(lpt_eles.pt),
                eta=ak.flatten(lpt_eles.eta),
                phi=ak.flatten(lpt_eles.phi))
    ele_kin.fill(sample=samp,
                ele_type="Candidate",
                pt=ak.flatten(cand_eles.pt),
                eta=ak.flatten(cand_eles.eta),
                phi=ak.flatten(cand_eles.phi))
    
    ele_trkHits.fill(sample=samp,
                    ele_type="Default",
                    numTrkHits=ak.flatten(eles.numTrackerHits),
                    numPixHits=ak.flatten(eles.numPixHits),
                    numStripHits=ak.flatten(eles.numStripHits))
    ele_trkHits.fill(sample=samp,
                    ele_type="Low pT",
                    numTrkHits=ak.flatten(lpt_eles.numTrackerHits),
                    numPixHits=ak.flatten(lpt_eles.numPixHits),
                    numStripHits=ak.flatten(lpt_eles.numStripHits))
    ele_trkHits.fill(sample=samp,
                    ele_type="Candidate",
                    numTrkHits=ak.flatten(cand_eles.numTrackerHits),
                    numPixHits=ak.flatten(cand_eles.numPixHits),
                    numStripHits=ak.flatten(cand_eles.numStripHits))

    ele_trkQual.fill(sample=samp,
                    ele_type="Default",
                    chi2=ak.flatten(eles.trkChi2),
                    trkIso=ak.flatten(eles.trkIso))
    ele_trkQual.fill(sample=samp,
                    ele_type="Low pT",
                    chi2=ak.flatten(lpt_eles.trkChi2),
                    trkIso=ak.flatten(lpt_eles.trkIso))
    ele_trkQual.fill(sample=samp,
                    ele_type="Candidate",
                    chi2=ak.flatten(cand_eles.trkChi2),
                    trkIso=ak.flatten(cand_eles.trkIso))
    
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
    gen_ele = events.GenEle
    gen_pos = events.GenPos

    ele_matches = events.GenEleMatches
    pos_matches = events.GenPosMatches
    ele_match = events.GenEleMatch
    pos_match = events.GenPosMatch
    ele_matchC = events.GenEleMatchC
    pos_matchC = events.GenPosMatchC

    # define labels for electron types & concatenate into larger collection
    ematch, pmatch = findUniqueMatches(ele_matches,ele_match,pos_matches,pos_match)
    ematch11, pmatch11 = findUniqueMatches(ele_matches,ele_match,pos_matches,pos_match,allow11=True)
    ematchC, pmatchC = findUniqueMatches(ele_matches,ele_matchC,pos_matches,pos_matchC,useCands=True)
    ematchC11, pmatchC11 = findUniqueMatches(ele_matches,ele_matchC,pos_matches,pos_matchC,allow11=True,useCands=True)

    # calculating other useful quantities
    nEmatch = ak.count(ele_matches[ele_matches.typ != 3].dr,axis=1)
    nPmatch = ak.count(pos_matches[pos_matches.typ != 3].dr,axis=1)
    nEmatchC = ak.count(ele_matches.dr,axis=1)
    nPmatchC = ak.count(pos_matches.dr,axis=1)

    # determining e+e- "match classes" for each event
    zeros = ak.zeros_like(ematch.typ) # generic vector of zeros for each event
    ones = ak.ones_like(ematch.typ) # generic vector fo ones for each event
    
    noMatch = (ematch11.typ == 0) & (pmatch11.typ == 0)
    oneMatch = ((ematch11.typ != 0) & (pmatch11.typ == 0)) | ((ematch11.typ == 0) & (pmatch11.typ != 0))
    degenMatch = (ematch11.typ != 0) & (pmatch11.typ != 0) & ((ematch11.typ == pmatch11.typ) & (ematch11.ind == pmatch11.ind))
    uniqMatch = (ematch11.typ != 0) & (pmatch11.typ != 0) & ((ematch11.typ != pmatch11.typ) | (ematch11.ind != pmatch11.ind))
    matchClass = ak.where(noMatch,zeros,-1*ones)
    matchClass = ak.where(oneMatch,ones,matchClass)
    matchClass = ak.where(degenMatch,2*ones,matchClass)
    matchClass = ak.where(uniqMatch,3*ones,matchClass)

    noMatchC = (ematchC11.typ == 0) & (pmatchC11.typ == 0)
    oneMatchC = ((ematchC11.typ != 0) & (pmatchC11.typ == 0)) | ((ematchC11.typ == 0) & (pmatchC11.typ != 0))
    degenMatchC = (ematchC11.typ != 0) & (pmatchC11.typ != 0) & ((ematchC11.typ == pmatchC11.typ) & (ematchC11.ind == pmatchC11.ind))
    uniqMatchC = (ematchC11.typ != 0) & (pmatchC11.typ != 0) & ((ematchC11.typ != pmatchC11.typ) | (ematchC11.ind != pmatchC11.ind))
    matchClassC = ak.where(noMatchC,zeros,-1*ones)
    matchClassC = ak.where(oneMatchC,ones,matchClassC)
    matchClassC = ak.where(degenMatchC,2*ones,matchClassC)
    matchClassC = ak.where(uniqMatchC,3*ones,matchClassC)

    # Fill match matrices, reg & low-pT matching
    match_mtx.fill(sample=samp,
                    set="RL",
                    scheme="unique",
                    Etype=ematch.typ,
                    Ptype=pmatch.typ)
    match_mtx.fill(sample=samp,
                    set="RL",
                    scheme="nearUnique",
                    Etype=ematch11.typ,
                    Ptype=pmatch11.typ)
    # Fill match matrices, reg & low-pT & cand matching
    match_mtx.fill(sample=samp,
                    set="RLC",
                    scheme="unique",
                    Etype=ematchC.typ,
                    Ptype=pmatchC.typ)
    match_mtx.fill(sample=samp,
                    set="RLC",
                    scheme="nearUnique",
                    Etype=ematchC11.typ,
                    Ptype=pmatchC11.typ)
    
    # Filling number-of-match matrices
    nMatch_mtx.fill(sample=samp,
                    set="RL",
                    nEmatch=nEmatch,
                    nPmatch=nPmatch)
    nMatch_mtx.fill(sample=samp,
                    set="RLC",
                    nEmatch=nEmatchC,
                    nPmatch=nPmatchC)

    # Filling match class matrices
    match_class.fill(sample=samp,
                    set="RL",
                    matchClass=matchClass)
    match_class.fill(sample=samp,
                    set="RLC",
                    matchClass=matchClassC)

    # Preparing variables and filling matched kinematic histograms
    mask_e1 = ematch.typ == 1
    mask_e2 = ematch.typ == 2
    mask_p1 = pmatch.typ == 1
    mask_p2 = pmatch.typ == 2
    matches_e1 = ak.flatten(events.Electron[mask_e1][ak.Array(ematch[mask_e1].ind[:,np.newaxis].to_list())])
    matches_e2 = ak.flatten(events.LptElectron[mask_e2][ak.Array(ematch[mask_e2].ind[:,np.newaxis].to_list())])
    matches_ele = ak.concatenate([matches_e1,matches_e2])
    matches_p1 = ak.flatten(events.Electron[mask_p1][ak.Array(pmatch[mask_p1].ind[:,np.newaxis].to_list())])
    matches_p2 = ak.flatten(events.LptElectron[mask_p2][ak.Array(pmatch[mask_p2].ind[:,np.newaxis].to_list())])
    matches_pos = ak.concatenate([matches_p1,matches_p2])
    matches = ak.concatenate([matches_ele,matches_pos])

    mask_e1C = ematchC.typ == 1
    mask_e2C = ematchC.typ == 2
    mask_e3C = ematchC.typ == 3
    mask_p1C = pmatchC.typ == 1
    mask_p2C = pmatchC.typ == 2
    mask_p3C = pmatchC.typ == 3
    matches_e1C = ak.flatten(events.Electron[mask_e1C][ak.Array(ematchC[mask_e1C].ind[:,np.newaxis].to_list())])
    matches_e2C = ak.flatten(events.LptElectron[mask_e2C][ak.Array(ematchC[mask_e2C].ind[:,np.newaxis].to_list())])
    matches_e3C = ak.flatten(events.EleCand[mask_e3C][ak.Array(ematchC[mask_e3C].ind[:,np.newaxis].to_list())])
    matches_eleC = ak.concatenate([matches_e1C,matches_e2C,matches_e3C])
    matches_p1C = ak.flatten(events.Electron[mask_p1C][ak.Array(pmatchC[mask_p1C].ind[:,np.newaxis].to_list())])
    matches_p2C = ak.flatten(events.LptElectron[mask_p2C][ak.Array(pmatchC[mask_p2C].ind[:,np.newaxis].to_list())])
    matches_p3C = ak.flatten(events.EleCand[mask_p3C][ak.Array(pmatchC[mask_p3C].ind[:,np.newaxis].to_list())])
    matches_posC = ak.concatenate([matches_p1C,matches_p2C,matches_p3C])
    matchesC = ak.concatenate([matches_eleC,matches_posC])

    matches_gen_ele = gen_ele[ematch.typ != 0]
    matches_gen_pos = gen_pos[pmatch.typ != 0]
    matches_gen = ak.concatenate([matches_gen_ele,matches_gen_pos])

    matches_gen_eleC = gen_ele[ematchC.typ != 0]
    matches_gen_posC = gen_pos[pmatchC.typ != 0]
    matches_genC = ak.concatenate([matches_gen_eleC,matches_gen_posC])

    # Filling with reg/lowpT matches
    match_ele_kin.fill(sample=samp,
                        set="RL",
                        pt=matches.pt,
                        eta=matches.eta,
                        phi=matches.phi)
    match_ele_gen_kin.fill(sample=samp,
                            set="RL",
                            pt=matches_gen.pt,
                            eta=matches_gen.eta,
                            phi=matches_gen.phi)
    match_ele_trkHits.fill(sample=samp,
                            set="RL",
                            numTrkHits=matches.numTrackerHits,
                            numPixHits=matches.numPixHits,
                            numStripHits=matches.numStripHits)
    match_ele_trkQual.fill(sample=samp,
                            set="RL",
                            chi2=matches.trkChi2,
                            trkIso=matches.trkIso)
    match_ele_gen_disp.fill(sample=samp,
                            set="RL",
                            vxy=matches_gen.vxy,
                            vz=matches_gen.vz)
    
    # Filling with reg/lowpT/cand matches
    match_ele_kin.fill(sample=samp,
                        set="RLC",
                        pt=matchesC.pt,
                        eta=matchesC.eta,
                        phi=matchesC.phi)
    match_ele_gen_kin.fill(sample=samp,
                            set="RLC",
                            pt=matches_genC.pt,
                            eta=matches_genC.eta,
                            phi=matches_genC.phi)
    match_ele_trkHits.fill(sample=samp,
                            set="RLC",
                            numTrkHits=matchesC.numTrackerHits,
                            numPixHits=matchesC.numPixHits,
                            numStripHits=matchesC.numStripHits)
    match_ele_trkQual.fill(sample=samp,
                            set="RLC",
                            chi2=matchesC.trkChi2,
                            trkIso=matchesC.trkIso)
    match_ele_gen_disp.fill(sample=samp,
                            set="RLC",
                            vxy=matches_genC.vxy,
                            vz=matches_genC.vz)
    
    return histos

def genParticles(events,histos,samp):
    gen_ele_disp = histos['gen_displacement']
    ele_kin = histos['ele_kinematics']
    ee_kin = histos["gen_ee_kinematics"]

    gen_eles = ak.concatenate([events.GenEle,events.GenPos])
    gen_ele = events.GenEle
    gen_pos = events.GenPos
    ee = events.genEE

    ele_kin.fill(sample=samp,ele_type="Generator",
                        pt=gen_eles.pt,eta=gen_eles.eta,phi=gen_eles.phi)
    gen_ele_disp.fill(sample=samp,
                        vxy=gen_eles.vxy,
                        vz=gen_eles.vz)
    ee_kin.fill(sample=samp,
                mass=ee.mass,
                dR=ee.dr)
    
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
    n_rr = ak.count(vtx_rr.vxy,axis=1)
    n_ll = ak.count(vtx_ll.vxy,axis=1)
    n_lr = ak.count(vtx_lr.vxy,axis=1)
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
    
    # Looking at displaced dilepton vertices found by standalone algorithms
    h_dEE_kin = histos["dispEE_vtx_kinematics"]
    h_dEE_qual = histos["dispEE_vtx_qual"]
    h_dEE_eff = histos["dispEE_eff"]

    dEE = events.EECand
    dEE_vxy = np.sqrt(dEE.vx**2 + dEE.vy**2)
    ndEE = ak.count(dEE.dR,axis=1)
    has_dEE = ndEE > 0

    # Masks for denominator categories
    noRegVtx = (n_rr == 0) & (n_lr == 0) & (n_ll == 0)

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

    return histos

def recoEE_vertices_genMatchByEles(events,histos,samp):
    # Reco vertices
    vtx_rr = events.RRvtx
    vtx_ll = events.LLvtx
    vtx_lr = events.LRvtx
    vtx_rc = events.RCvtx
    vtx_lc = events.LCvtx
    vtx_cc = events.dispEE

    ele = events.Electron
    lpt_ele = events.LptElectron
    cand_ele = events.EleCand
    gen_ele = events.GenEle
    gen_pos = events.GenPos
    ele_matches = events.GenEleMatches
    ele_match = events.GenEleMatch
    ele_matchC = events.GenEleMatchC
    pos_matches = events.GenPosMatches
    pos_match = events.GenPosMatch
    pos_matchC = events.GenPosMatchC

    ematch, pmatch = findUniqueMatches(ele_matches,ele_match,pos_matches,pos_match)
    ematchC, pmatchC = findUniqueMatches(ele_matches,ele_matchC,pos_matches,pos_matchC,useCands=True)

    # masks for different vertex types, low/reg only
    noMatch_RL = (ematch.typ == 0) & (pmatch.typ == 0)
    matchRR_RL = (ematch.typ == 1) & (pmatch.typ == 1)
    matchLR_RL = ((ematch.typ == 1) & (pmatch.typ == 2)) | ((ematch.typ == 2) & (pmatch.typ == 1))
    matchLL_RL = (ematch.typ == 2) & (pmatch.typ == 2)
    ones = ak.ones_like(ematch.typ)
    zeros = ak.zeros_like(ematch.typ)
    matchType_RL = ak.where(noMatch_RL,zeros,-1*ones)
    matchType_RL = ak.where(matchRR_RL,ones,matchType_RL)
    matchType_RL = ak.where(matchLR_RL,2*ones,matchType_RL)
    matchType_RL = ak.where(matchLL_RL,3*ones,matchType_RL)

    # masks for different vertex types, low/reg/cand
    noMatch_RLC = (ematchC.typ == 0) & (pmatchC.typ == 0)
    matchRR_RLC = (ematchC.typ == 1) & (pmatchC.typ == 1)
    matchRC_RLC = ((ematchC.typ == 1) & (pmatchC.typ == 3)) | ((ematchC.typ == 3) & (pmatchC.typ == 1))
    matchLR_RLC = ((ematchC.typ == 1) & (pmatchC.typ == 2)) | ((ematchC.typ == 2) & (pmatchC.typ == 1))
    matchLC_RLC = ((ematchC.typ == 2) & (pmatchC.typ == 3)) | ((ematchC.typ == 3) & (pmatchC.typ == 2))
    matchLL_RLC = (ematchC.typ == 2) & (pmatchC.typ == 2)
    matchCC_RLC = (ematchC.typ == 3) & (pmatchC.typ == 3)
    ones = ak.ones_like(ematchC.typ)
    zeros = ak.zeros_like(ematchC.typ)
    matchType_RLC = ak.where(noMatch_RLC,zeros,-1*ones)
    matchType_RLC = ak.where(matchRR_RLC,ones,matchType_RLC)
    matchType_RLC = ak.where(matchLR_RLC,2*ones,matchType_RLC)
    matchType_RLC = ak.where(matchLL_RLC,3*ones,matchType_RLC)
    matchType_RLC = ak.where(matchRC_RLC,4*ones,matchType_RLC)
    matchType_RLC = ak.where(matchLC_RLC,5*ones,matchType_RLC)
    matchType_RLC = ak.where(matchCC_RLC,6*ones,matchType_RLC)

    # reg/lowpT only
    RRvtx_RL = vtx_rr[matchRR_RL]
    RRvtx_mask_RL = ((ematch[matchRR_RL].ind == RRvtx_RL.idx1) & (pmatch[matchRR_RL].ind == RRvtx_RL.idx2)) |\
                    ((ematch[matchRR_RL].ind == RRvtx_RL.idx2) & (pmatch[matchRR_RL].ind == RRvtx_RL.idx1))
    RRvtx_RL = ak.flatten(RRvtx_RL[RRvtx_mask_RL])

    LRvtx_RL = vtx_lr[matchLR_RL]
    LRvtx_mask_RL = ((ematch[matchLR_RL].typ == 2) & (ematch[matchLR_RL].ind == LRvtx_RL.idx1) & (pmatch[matchLR_RL].ind == LRvtx_RL.idx2)) |\
                    ((ematch[matchLR_RL].typ == 1) & (ematch[matchLR_RL].ind == LRvtx_RL.idx2) & (pmatch[matchLR_RL].ind == LRvtx_RL.idx1))
    LRvtx_RL = ak.flatten(LRvtx_RL[LRvtx_mask_RL])

    LLvtx_RL = vtx_ll[matchLL_RL]
    LLvtx_mask_RL = ((ematch[matchLL_RL].ind == LLvtx_RL.idx1) & (pmatch[matchLL_RL].ind == LLvtx_RL.idx2)) |\
                    ((ematch[matchLL_RL].ind == LLvtx_RL.idx2) & (pmatch[matchLL_RL].ind == LLvtx_RL.idx1))
    LLvtx_RL = ak.flatten(LLvtx_RL[LLvtx_mask_RL])

    # reg,lowpT,cand
    RRvtx_RLC = vtx_rr[matchRR_RLC]
    RRvtx_mask_RLC = ((ematchC[matchRR_RLC].ind == RRvtx_RLC.idx1) & (pmatchC[matchRR_RLC].ind == RRvtx_RLC.idx2)) |\
                     ((ematchC[matchRR_RLC].ind == RRvtx_RLC.idx2) & (pmatchC[matchRR_RLC].ind == RRvtx_RLC.idx1))
    RRvtx_RLC = ak.flatten(RRvtx_RLC[RRvtx_mask_RLC])

    RCvtx_RLC = vtx_rc[matchRC_RLC]
    RCvtx_mask_RLC = ((ematchC[matchRC_RLC].typ == 1) & (ematchC[matchRC_RLC].ind == RCvtx_RLC.idx1) & (pmatchC[matchRC_RLC].ind == RCvtx_RLC.idx2)) |\
                     ((ematchC[matchRC_RLC].typ == 3) & (ematchC[matchRC_RLC].ind == RCvtx_RLC.idx2) & (pmatchC[matchRC_RLC].ind == RCvtx_RLC.idx1))
    RCvtx_RLC = ak.flatten(RCvtx_RLC[RCvtx_mask_RLC])

    LCvtx_RLC = vtx_lc[matchLC_RLC]
    LCvtx_mask_RLC = ((ematchC[matchLC_RLC].typ == 2) & (ematchC[matchLC_RLC].ind == LCvtx_RLC.idx1) & (pmatchC[matchLC_RLC].ind == LCvtx_RLC.idx2)) |\
                     ((ematchC[matchLC_RLC].typ == 3) & (ematchC[matchLC_RLC].ind == LCvtx_RLC.idx2) & (pmatchC[matchLC_RLC].ind == LCvtx_RLC.idx1))
    LCvtx_RLC = ak.flatten(LCvtx_RLC[LCvtx_mask_RLC])

    LRvtx_RLC = vtx_lr[matchLR_RLC]
    LRvtx_mask_RLC = ((ematchC[matchLR_RLC].typ == 2) & (ematchC[matchLR_RLC].ind == LRvtx_RLC.idx1) & (pmatchC[matchLR_RLC].ind == LRvtx_RLC.idx2)) |\
                     ((ematchC[matchLR_RLC].typ == 1) & (ematchC[matchLR_RLC].ind == LRvtx_RLC.idx2) & (pmatchC[matchLR_RLC].ind == LRvtx_RLC.idx1))
    LRvtx_RLC = ak.flatten(LRvtx_RLC[LRvtx_mask_RLC])

    LLvtx_RLC = vtx_ll[matchLL_RLC]
    LLvtx_mask_RLC = ((ematchC[matchLL_RLC].ind == LLvtx_RLC.idx1) & (pmatchC[matchLL_RLC].ind == LLvtx_RLC.idx2)) |\
                     ((ematchC[matchLL_RLC].ind == LLvtx_RLC.idx2) & (pmatchC[matchLL_RLC].ind == LLvtx_RLC.idx1))
    LLvtx_RLC = ak.flatten(LLvtx_RLC[LLvtx_mask_RLC])

    # Filling histograms, finally!
    h_vtx = histos["vtx_genMatchByEle"]
    h_vtx_kin = histos["vtx_kin_genMatchByEle"]
    h_vtx_disp = histos["vtx_disp_genMatchByEle"]
    h_vtx_type = histos["vtx_matchTypeByEle"]

    # filling with the reco/lowpT set first
    h_vtx_type.fill(sample=samp,
                    set="RL",
                    type=matchType_RL)

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
    h_vtx_type.fill(sample=samp,
                    set="RLC",
                    type=matchType_RLC)

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

    # Filling EE candidate info
    # Events with no gen-matched electrons (reg or low-pt)
    h_dEE_kin = histos["dispEE_vtx_kinematics"]
    h_dEE_qual = histos["dispEE_vtx_qual"]
    h_dEE_eff = histos["dispEE_eff"]

    dEE = events.EECand
    dEE_vxy = np.sqrt(dEE.vx**2 + dEE.vy**2)
    ndEE = ak.count(dEE.dR,axis=1)
    has_dEE = ndEE > 0
    noGenMatchVtx = (ematch.typ == 0) | (pmatch.typ == 0)


    h_dEE_eff.fill(sample=samp,
                    denom="noGenMatchVtx",
                    ndEE=ndEE[noGenMatchVtx])
    h_dEE_kin.fill(sample=samp,
                    denom="noGenMatchVtx",
                    mass=ak.flatten(dEE[has_dEE & noGenMatchVtx].mass),
                    leadPt=ak.flatten(dEE[has_dEE & noGenMatchVtx].leadingPt),
                    dR=ak.flatten(dEE[has_dEE & noGenMatchVtx].dR))
    h_dEE_qual.fill(sample=samp,
                    denom="noGenMatchVtx",
                    vxy=ak.flatten(dEE_vxy[has_dEE & noGenMatchVtx]),
                    dxy=ak.flatten(dEE[has_dEE & noGenMatchVtx].trackDxy_PV),
                    chi2=ak.flatten(dEE[has_dEE & noGenMatchVtx].normalizedChi2))
    
    return histos
from histobins import *
from coffea.hist import Hist
import numpy as np
import awkward as ak

histograms = {
    "ele_kinematics" : Hist("Events",sample,cut,ele_type,ele_pt,ele_eta,ele_phi),
    "gen_ee_kinematics" : Hist("Events",sample,cut,mass,dR),
    "ele_match_matrix" : Hist("Events",sample,cut,ele_set,scheme,etype,ptype),
    "ele_nMatch_matrix" : Hist("Events",sample,cut,ele_set,nEmatch,nPmatch),
    "ele_match_class" : Hist("Events",sample,cut,ele_set,matchClass),
    "matched_ele_kinematics" : Hist("Events",sample,cut,ele_set,ele_pt,ele_eta,ele_phi),
    "matched_ele_gen_kinematics" : Hist("Events",sample,cut,ele_set,ele_pt,ele_eta,ele_phi),
    "ele_trackHits" : Hist("Events",sample,cut,ele_type,ele_trkHits,ele_pixHits,ele_stripHits),
    "matched_ele_trackHits" : Hist("Events",sample,cut,ele_set,ele_trkHits,ele_pixHits,ele_stripHits),
    "ele_trackQual" : Hist("Events",sample,cut,ele_type,ele_chi2,ele_trkIso),
    "matched_ele_trackQual" : Hist("Events",sample,cut,ele_set,ele_chi2,ele_trkIso),
    "matched_ele_trackRes" : Hist("Events",sample,cut,ele_set,ele_angRes),
    "gen_displacement" : Hist("Events",sample,cut,vxy,vz),
    "matched_ele_gen_displacement" : Hist("Events",sample,cut,ele_set,vxy,vz),
    "RRvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "RRvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
    "LLvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "LLvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
    "LRvtx_vxy" : Hist("Events",sample,cut,vxy,vtx_vxyErr),
    "LRvtx_stats" : Hist("Events",sample,cut,vtx_chi2,vtx_vxySignif),
    "vtx_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,vtx_chi2,dR),
    "vtx_qual_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,vtx_chi2,vtx_prob),
    "vtx_chi2_vs_vxy_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,vtx_chi2,vxy),
    "vtx_prob_vs_vxy_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,vxy,vtx_prob),
    "vtx_disp_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,vxy,vtx_vxySignif),
    "vtx_kin_genMatchByEle" : Hist("Events",sample,cut,vtx_type,ele_set,mass,energy),
    "vtx_matchTypeByEle" : Hist("Events",sample,cut,ele_set,vtx_matchType),
    "match_matrix" : Hist("Events",sample,cut,etype,ptype)
}

subroutines = ["eleVtxGenMatching"]

def fillHistos(events,histos,samp,cut):
    # Electron information
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Default",
                pt=ak.flatten(events.Electron.pt),
                eta=ak.flatten(events.Electron.eta),
                phi=ak.flatten(events.Electron.phi))
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Low pT",
                pt=ak.flatten(events.LptElectron.pt),
                eta=ak.flatten(events.LptElectron.eta),
                phi=ak.flatten(events.LptElectron.phi))
    histos["ele_kinematics"].fill(sample=samp,cut=cut,
                ele_type="Candidate",
                pt=ak.flatten(events.EleCand.pt),
                eta=ak.flatten(events.EleCand.eta),
                phi=ak.flatten(events.EleCand.phi))
    
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Default",
                    numTrkHits=ak.flatten(events.Electron.numTrackerHits),
                    numPixHits=ak.flatten(events.Electron.numPixHits),
                    numStripHits=ak.flatten(events.Electron.numStripHits))
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Low pT",
                    numTrkHits=ak.flatten(events.LptElectron.numTrackerHits),
                    numPixHits=ak.flatten(events.LptElectron.numPixHits),
                    numStripHits=ak.flatten(events.LptElectron.numStripHits))
    histos["ele_trackHits"].fill(sample=samp,cut=cut,
                    ele_type="Candidate",
                    numTrkHits=ak.flatten(events.EleCand.numTrackerHits),
                    numPixHits=ak.flatten(events.EleCand.numPixHits),
                    numStripHits=ak.flatten(events.EleCand.numStripHits))

    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Default",
                    chi2=ak.flatten(events.Electron.trkChi2),
                    trkIso=ak.flatten(events.Electron.trkIso))
    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Low pT",
                    chi2=ak.flatten(events.LptElectron.trkChi2),
                    trkIso=ak.flatten(events.LptElectron.trkIso))
    histos['ele_trackQual'].fill(sample=samp,cut=cut,
                    ele_type="Candidate",
                    chi2=ak.flatten(events.EleCand.trkChi2),
                    trkIso=ak.flatten(events.EleCand.trkIso))

    # Fill match matrices, reg & low-pT matching
    histos["ele_match_matrix"].fill(sample=samp,cut=cut,
                    set="RL",
                    scheme="unique",
                    Etype=events.ematch.typ,
                    Ptype=events.pmatch.typ)
    histos["ele_match_matrix"].fill(sample=samp,cut=cut,
                    set="RL",
                    scheme="nearUnique",
                    Etype=events.ematch11.typ,
                    Ptype=events.pmatch11.typ)
    # Fill match matrices, reg & low-pT & cand matching
    histos["ele_match_matrix"].fill(sample=samp,cut=cut,
                    set="RLC",
                    scheme="unique",
                    Etype=events.ematchC.typ,
                    Ptype=events.pmatchC.typ)
    histos["ele_match_matrix"].fill(sample=samp,cut=cut,
                    set="RLC",
                    scheme="nearUnique",
                    Etype=events.ematchC11.typ,
                    Ptype=events.pmatchC11.typ)
    
    # Filling number-of-match matrices
    histos["ele_nMatch_matrix"].fill(sample=samp,cut=cut,
                    set="RL",
                    nEmatch=events.nEmatch,
                    nPmatch=events.nPmatch)
    histos["ele_nMatch_matrix"].fill(sample=samp,cut=cut,
                    set="RLC",
                    nEmatch=events.nEmatchC,
                    nPmatch=events.nPmatchC)

    # Filling match class matrices
    histos["ele_match_class"].fill(sample=samp,cut=cut,
                    set="RL",
                    matchClass=events.eleMatchClass)
    histos["ele_match_class"].fill(sample=samp,cut=cut,
                    set="RLC",
                    matchClass=events.eleMatchClassC)

    # Filling with reg/lowpT matches
    mask_e1 = events.ematch.typ == 1
    mask_e2 = events.ematch.typ == 2
    mask_p1 = events.pmatch.typ == 1
    mask_p2 = events.pmatch.typ == 2
    matches_e1 = ak.flatten(events.Electron[mask_e1][ak.Array(events.ematch[mask_e1].ind[:,np.newaxis].to_list())])
    matches_e2 = ak.flatten(events.LptElectron[mask_e2][ak.Array(events.ematch[mask_e2].ind[:,np.newaxis].to_list())])
    matches_ele = ak.concatenate([matches_e1,matches_e2])
    matches_p1 = ak.flatten(events.Electron[mask_p1][ak.Array(events.pmatch[mask_p1].ind[:,np.newaxis].to_list())])
    matches_p2 = ak.flatten(events.LptElectron[mask_p2][ak.Array(events.pmatch[mask_p2].ind[:,np.newaxis].to_list())])
    matches_pos = ak.concatenate([matches_p1,matches_p2])
    matches = ak.concatenate([matches_ele,matches_pos])

    mask_e1C = events.ematchC.typ == 1
    mask_e2C = events.ematchC.typ == 2
    mask_e3C = events.ematchC.typ == 3
    mask_p1C = events.pmatchC.typ == 1
    mask_p2C = events.pmatchC.typ == 2
    mask_p3C = events.pmatchC.typ == 3
    matches_e1C = ak.flatten(events.Electron[mask_e1C][ak.Array(events.ematchC[mask_e1C].ind[:,np.newaxis].to_list())])
    matches_e2C = ak.flatten(events.LptElectron[mask_e2C][ak.Array(events.ematchC[mask_e2C].ind[:,np.newaxis].to_list())])
    matches_e3C = ak.flatten(events.EleCand[mask_e3C][ak.Array(events.ematchC[mask_e3C].ind[:,np.newaxis].to_list())])
    matches_eleC = ak.concatenate([matches_e1C,matches_e2C,matches_e3C])
    matches_p1C = ak.flatten(events.Electron[mask_p1C][ak.Array(events.pmatchC[mask_p1C].ind[:,np.newaxis].to_list())])
    matches_p2C = ak.flatten(events.LptElectron[mask_p2C][ak.Array(events.pmatchC[mask_p2C].ind[:,np.newaxis].to_list())])
    matches_p3C = ak.flatten(events.EleCand[mask_p3C][ak.Array(events.pmatchC[mask_p3C].ind[:,np.newaxis].to_list())])
    matches_posC = ak.concatenate([matches_p1C,matches_p2C,matches_p3C])
    matchesC = ak.concatenate([matches_eleC,matches_posC])

    matches_gen_ele = events.GenEle[events.ematch.typ != 0]
    matches_gen_pos = events.GenPos[events.pmatch.typ != 0]
    matches_gen = ak.concatenate([matches_gen_ele,matches_gen_pos])

    matches_gen_eleC = events.GenEle[events.ematchC.typ != 0]
    matches_gen_posC = events.GenPos[events.pmatchC.typ != 0]
    matches_genC = ak.concatenate([matches_gen_eleC,matches_gen_posC])

    histos['matched_ele_kinematics'].fill(sample=samp,cut=cut,
                        set="RL",
                        pt=matches.pt,
                        eta=matches.eta,
                        phi=matches.phi)
    histos['matched_ele_gen_kinematics'].fill(sample=samp,cut=cut,
                            set="RL",
                            pt=matches_gen.pt,
                            eta=matches_gen.eta,
                            phi=matches_gen.phi)
    histos['matched_ele_trackHits'].fill(sample=samp,cut=cut,
                            set="RL",
                            numTrkHits=matches.numTrackerHits,
                            numPixHits=matches.numPixHits,
                            numStripHits=matches.numStripHits)
    histos['matched_ele_trackQual'].fill(sample=samp,cut=cut,
                            set="RL",
                            chi2=matches.trkChi2,
                            trkIso=matches.trkIso)
    histos['matched_ele_gen_displacement'].fill(sample=samp,cut=cut,
                            set="RL",
                            vxy=matches_gen.vxy,
                            vz=matches_gen.vz)
    
    # Filling with reg/lowpT/cand matches
    histos['matched_ele_kinematics'].fill(sample=samp,cut=cut,
                        set="RLC",
                        pt=matchesC.pt,
                        eta=matchesC.eta,
                        phi=matchesC.phi)
    histos['matched_ele_gen_kinematics'].fill(sample=samp,cut=cut,
                            set="RLC",
                            pt=matches_genC.pt,
                            eta=matches_genC.eta,
                            phi=matches_genC.phi)
    histos['matched_ele_trackHits'].fill(sample=samp,cut=cut,
                            set="RLC",
                            numTrkHits=matchesC.numTrackerHits,
                            numPixHits=matchesC.numPixHits,
                            numStripHits=matchesC.numStripHits)
    histos['matched_ele_trackQual'].fill(sample=samp,cut=cut,
                            set="RLC",
                            chi2=matchesC.trkChi2,
                            trkIso=matchesC.trkIso)
    histos['matched_ele_gen_displacement'].fill(sample=samp,cut=cut,
                            set="RLC",
                            vxy=matches_genC.vxy,
                            vz=matches_genC.vz)

    # Filling gen particle information (signal only)
    gen_eles = ak.concatenate([events.GenEle,events.GenPos])
    histos['ele_kinematics'].fill(sample=samp,cut=cut,ele_type="Generator",
                        pt=gen_eles.pt,eta=gen_eles.eta,phi=gen_eles.phi)
    histos['gen_displacement'].fill(sample=samp,cut=cut,
                        vxy=gen_eles.vxy,
                        vz=gen_eles.vz)
    histos["gen_ee_kinematics"].fill(sample=samp,cut=cut,
                mass=events.genEE.mass,
                dR=events.genEE.dr)

    # Reco vertex information
    histos["RRvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.RRvtx.vxy),
                        sigma_vxy=ak.flatten(events.RRvtx.sigmavxy))
    histos["RRvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.RRvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.RRvtx.vxy/events.RRvtx.sigmavxy))
    # low-low, i.e. two low-pT electrons matched to a vertex
    histos["LLvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.LLvtx.vxy),
                        sigma_vxy=ak.flatten(events.LLvtx.sigmavxy))
    histos["LLvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.LLvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.LLvtx.vxy/events.LLvtx.sigmavxy))
    # low-reg, i.e. one low-pT and one default electron matched to a vertex
    histos["LRvtx_vxy"].fill(sample=samp,cut=cut,
                        vxy=ak.flatten(events.LRvtx.vxy),
                        sigma_vxy=ak.flatten(events.LRvtx.sigmavxy))
    histos["LRvtx_stats"].fill(sample=samp,cut=cut,
                                chi2=ak.flatten(events.LRvtx.reduced_chi2),
                            vxy_signif=ak.flatten(events.LRvtx.vxy/events.LRvtx.sigmavxy))

    # Gen-matched vertices from gen-matched electrons
    RRvtx_RL = events.RRvtx[events.eeMatchRR_RL]
    RRvtx_mask_RL = ((events.ematch[events.eeMatchRR_RL].ind == RRvtx_RL.idx1) & (events.pmatch[events.eeMatchRR_RL].ind == RRvtx_RL.idx2)) |\
                    ((events.ematch[events.eeMatchRR_RL].ind == RRvtx_RL.idx2) & (events.pmatch[events.eeMatchRR_RL].ind == RRvtx_RL.idx1))
    RRvtx_RL = ak.flatten(RRvtx_RL[RRvtx_mask_RL])

    LRvtx_RL = events.LRvtx[events.eeMatchLR_RL]
    LRvtx_mask_RL = ((events.ematch[events.eeMatchLR_RL].typ == 2) & (events.ematch[events.eeMatchLR_RL].ind == LRvtx_RL.idx1) & (events.pmatch[events.eeMatchLR_RL].ind == LRvtx_RL.idx2)) |\
                    ((events.ematch[events.eeMatchLR_RL].typ == 1) & (events.ematch[events.eeMatchLR_RL].ind == LRvtx_RL.idx2) & (events.pmatch[events.eeMatchLR_RL].ind == LRvtx_RL.idx1))
    LRvtx_RL = ak.flatten(LRvtx_RL[LRvtx_mask_RL])

    LLvtx_RL = events.LLvtx[events.eeMatchLL_RL]
    LLvtx_mask_RL = ((events.ematch[events.eeMatchLL_RL].ind == LLvtx_RL.idx1) & (events.pmatch[events.eeMatchLL_RL].ind == LLvtx_RL.idx2)) |\
                    ((events.ematch[events.eeMatchLL_RL].ind == LLvtx_RL.idx2) & (events.pmatch[events.eeMatchLL_RL].ind == LLvtx_RL.idx1))
    LLvtx_RL = ak.flatten(LLvtx_RL[LLvtx_mask_RL])

    RRvtx_RLC = events.RRvtx[events.eeMatchRR_RLC]
    RRvtx_mask_RLC = ((events.ematchC[events.eeMatchRR_RLC].ind == RRvtx_RLC.idx1) & (events.pmatchC[events.eeMatchRR_RLC].ind == RRvtx_RLC.idx2)) |\
                     ((events.ematchC[events.eeMatchRR_RLC].ind == RRvtx_RLC.idx2) & (events.pmatchC[events.eeMatchRR_RLC].ind == RRvtx_RLC.idx1))
    RRvtx_RLC = ak.flatten(RRvtx_RLC[RRvtx_mask_RLC])

    RCvtx_RLC = events.RCvtx[events.eeMatchRC_RLC]
    RCvtx_mask_RLC = ((events.ematchC[events.eeMatchRC_RLC].typ == 1) & (events.ematchC[events.eeMatchRC_RLC].ind == RCvtx_RLC.idx1) & (events.pmatchC[events.eeMatchRC_RLC].ind == RCvtx_RLC.idx2)) |\
                     ((events.ematchC[events.eeMatchRC_RLC].typ == 3) & (events.ematchC[events.eeMatchRC_RLC].ind == RCvtx_RLC.idx2) & (events.pmatchC[events.eeMatchRC_RLC].ind == RCvtx_RLC.idx1))
    RCvtx_RLC = ak.flatten(RCvtx_RLC[RCvtx_mask_RLC])

    LCvtx_RLC = events.LCvtx[events.eeMatchLC_RLC]
    LCvtx_mask_RLC = ((events.ematchC[events.eeMatchLC_RLC].typ == 2) & (events.ematchC[events.eeMatchLC_RLC].ind == LCvtx_RLC.idx1) & (events.pmatchC[events.eeMatchLC_RLC].ind == LCvtx_RLC.idx2)) |\
                     ((events.ematchC[events.eeMatchLC_RLC].typ == 3) & (events.ematchC[events.eeMatchLC_RLC].ind == LCvtx_RLC.idx2) & (events.pmatchC[events.eeMatchLC_RLC].ind == LCvtx_RLC.idx1))
    LCvtx_RLC = ak.flatten(LCvtx_RLC[LCvtx_mask_RLC])

    LRvtx_RLC = events.LRvtx[events.eeMatchLR_RLC]
    LRvtx_mask_RLC = ((events.ematchC[events.eeMatchLR_RLC].typ == 2) & (events.ematchC[events.eeMatchLR_RLC].ind == LRvtx_RLC.idx1) & (events.pmatchC[events.eeMatchLR_RLC].ind == LRvtx_RLC.idx2)) |\
                     ((events.ematchC[events.eeMatchLR_RLC].typ == 1) & (events.ematchC[events.eeMatchLR_RLC].ind == LRvtx_RLC.idx2) & (events.pmatchC[events.eeMatchLR_RLC].ind == LRvtx_RLC.idx1))
    LRvtx_RLC = ak.flatten(LRvtx_RLC[LRvtx_mask_RLC])

    LLvtx_RLC = events.LLvtx[events.eeMatchLL_RLC]
    LLvtx_mask_RLC = ((events.ematchC[events.eeMatchLL_RLC].ind == LLvtx_RLC.idx1) & (events.pmatchC[events.eeMatchLL_RLC].ind == LLvtx_RLC.idx2)) |\
                     ((events.ematchC[events.eeMatchLL_RLC].ind == LLvtx_RLC.idx2) & (events.pmatchC[events.eeMatchLL_RLC].ind == LLvtx_RLC.idx1))
    LLvtx_RLC = ak.flatten(LLvtx_RLC[LLvtx_mask_RLC])

    histos["vtx_matchTypeByEle"].fill(sample=samp,cut=cut,
                    set="RL",
                    type=events.eeMatchType_RL)

    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Reg-Reg",
               set="RL",
               chi2=RRvtx_RL.reduced_chi2,
               dR=RRvtx_RL.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Low-Reg",
               set="RL",
               chi2=LRvtx_RL.reduced_chi2,
               dR=LRvtx_RL.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Low-Low",
               set="RL",
               chi2=LLvtx_RL.reduced_chi2,
               dR=LLvtx_RL.dR)
    
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Reg-Reg",
                   set="RL",
                   mass=RRvtx_RL.m,
                   energy=RRvtx_RL.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Low-Reg",
                   set="RL",
                   mass=LRvtx_RL.m,
                   energy=LRvtx_RL.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Low-Low",
                   set="RL",
                   mass=LLvtx_RL.m,
                   energy=LLvtx_RL.energy)
    
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Reg",
                    set="RL",
                    vxy=RRvtx_RL.vxy,
                    vxy_signif=(RRvtx_RL.vxy/RRvtx_RL.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Reg",
                    set="RL",
                    vxy=LRvtx_RL.vxy,
                    vxy_signif=(LRvtx_RL.vxy/LRvtx_RL.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Low",
                    set="RL",
                    vxy=LLvtx_RL.vxy,
                    vxy_signif=(LLvtx_RL.vxy/LLvtx_RL.sigmavxy))

    histos["vtx_qual_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Reg",
                    set="RL",
                    prob=RRvtx_RL.prob,
                    chi2=RRvtx_RL.reduced_chi2)
    histos["vtx_qual_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Reg",
                    set="RL",
                    prob=LRvtx_RL.prob,
                    chi2=LRvtx_RL.reduced_chi2)
    histos["vtx_qual_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Low",
                    set="RL",
                    prob=LLvtx_RL.prob,
                    chi2=LLvtx_RL.reduced_chi2)

    histos["vtx_chi2_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Reg",
                    set="RL",
                    vxy=RRvtx_RL.vxy,
                    chi2=RRvtx_RL.reduced_chi2)
    histos["vtx_chi2_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Reg",
                    set="RL",
                    vxy=LRvtx_RL.vxy,
                    chi2=LRvtx_RL.reduced_chi2)
    histos["vtx_chi2_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Low",
                    set="RL",
                    vxy=LLvtx_RL.vxy,
                    chi2=LLvtx_RL.reduced_chi2)

    histos["vtx_prob_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Reg",
                    set="RL",
                    vxy=RRvtx_RL.vxy,
                    prob=RRvtx_RL.prob)
    histos["vtx_prob_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Reg",
                    set="RL",
                    vxy=LRvtx_RL.vxy,
                    prob=LRvtx_RL.prob)
    histos["vtx_prob_vs_vxy_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Low",
                    set="RL",
                    vxy=LLvtx_RL.vxy,
                    prob=LLvtx_RL.prob)

    # Filling with reg/lowpT/cand set
    histos["vtx_matchTypeByEle"].fill(sample=samp,cut=cut,
                    set="RLC",
                    type=events.eeMatchType_RLC)

    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Reg-Reg",
               set="RLC",
               chi2=RRvtx_RLC.reduced_chi2,
               dR=RRvtx_RLC.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Reg-Cand",
               set="RLC",
               chi2=RCvtx_RLC.reduced_chi2,
               dR=RCvtx_RLC.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Low-Reg",
               set="RLC",
               chi2=LRvtx_RLC.reduced_chi2,
               dR=LRvtx_RLC.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Low-Cand",
               set="RLC",
               chi2=LCvtx_RLC.reduced_chi2,
               dR=LCvtx_RLC.dR)
    histos["vtx_genMatchByEle"].fill(sample=samp,cut=cut,
               type="Low-Low",
               set="RLC",
               chi2=LLvtx_RLC.reduced_chi2,
               dR=LLvtx_RLC.dR)
    
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Reg-Reg",
                   set="RLC",
                   mass=RRvtx_RLC.m,
                   energy=RRvtx_RLC.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Reg-Cand",
                   set="RLC",
                   mass=RCvtx_RLC.m,
                   energy=RCvtx_RLC.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Low-Reg",
                   set="RLC",
                   mass=LRvtx_RLC.m,
                   energy=LRvtx_RLC.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Low-Cand",
                   set="RLC",
                   mass=LCvtx_RLC.m,
                   energy=LCvtx_RLC.energy)
    histos["vtx_kin_genMatchByEle"].fill(sample=samp,cut=cut,
                   type="Low-Low",
                   set="RLC",
                   mass=LLvtx_RLC.m,
                   energy=LLvtx_RLC.energy)
    
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Reg",
                    set="RLC",
                    vxy=RRvtx_RLC.vxy,
                    vxy_signif=(RRvtx_RLC.vxy/RRvtx_RLC.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Reg-Cand",
                    set="RLC",
                    vxy=RCvtx_RLC.vxy,
                    vxy_signif=(RCvtx_RLC.vxy/RCvtx_RLC.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Reg",
                    set="RLC",
                    vxy=LRvtx_RLC.vxy,
                    vxy_signif=(LRvtx_RLC.vxy/LRvtx_RLC.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Cand",
                    set="RLC",
                    vxy=LCvtx_RLC.vxy,
                    vxy_signif=(LCvtx_RLC.vxy/LCvtx_RLC.sigmavxy))
    histos["vtx_disp_genMatchByEle"].fill(sample=samp,cut=cut,
                    type="Low-Low",
                    set="RLC",
                    vxy=LLvtx_RLC.vxy,
                    vxy_signif=(LLvtx_RLC.vxy/LLvtx_RLC.sigmavxy))
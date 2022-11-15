import awkward as ak
import numpy as np

def mergeTrees(files):
    ## samples = ["Mchi-42p0_dMchi-4p0_ctau-10","Mchi-48p0_dMchi-16p0_ctau-1","Mchi-55p0_dMchi-10p0_ctau-1"]
    e1 = []
    le1 = []
    genE1 = []
    genP1 = []
    genEMatch1 = []
    genPMatch1 = []
    genEMatches1 = []
    genPMatches1 = []
    for fname in files:
        t1 = tools.loadNano(fname)
        n1 = t1.eventNum
        o1 = ak.argsort(n1)
        e1.append(t1.Electron[o1])
        le1.append(t1.LptElectron[o1])
        genE1.append(t1.GenEle[o1])
        genP1.append(t1.GenPos[o1])
        genEMatch1.append(t1.GenEleMatch[o1])
        genPMatch1.append(t1.GenPosMatch[o1])
        genEMatches1.append(t1.GenEleMatches[o1])
        genPMatches1.append(t1.GenPosMatches[o1])

    e1 = ak.concatenate(e1)
    le1 = ak.concatenate(le1)
    genE1 = ak.concatenate(genE1)
    genP1 = ak.concatenate(genP1)
    genEMatch1 = ak.concatenate(genEMatch1)
    genPMatch1 = ak.concatenate(genPMatch1)
    genEMatches1 = ak.concatenate(genEMatches1)
    genPMatches1 = ak.concatenate(genPMatches1)

    return e1,le1,genE1,genP1,genEMatch1,genPMatch1,genEMatches1,genpMatches1

def getMatchesv2(tree,eType,dRcut=0.3):
    ie = tree.GenEleClosest.ind
    te = tree.GenEleClosest.typ
    re = tree.GenEleClosest.dr

    ip = tree.GenPosClosest.ind
    tp = tree.GenPosClosest.typ
    rp = tree.GenPosClosest.dr

    eMatch = (te==eType) & (re<dRcut)
    ieMatch = ak.Array(ie[eMatch][:,np.newaxis].to_list())

    pMatch = (tp==eType) & (rp<dRcut)
    ipMatch = ak.Array(ip[pMatch][:,np.newaxis].to_list())

    electrons = tree.LptElectron if eType == 2 else tree.Electron
    
    genEmatch = tree.GenEle[eMatch]
    recoEmatch = electrons[eMatch][ieMatch][:,0]

    genPmatch = tree.GenPos[pMatch]
    recoPmatch = electrons[pMatch][ipMatch][:,0]
    
    return ak.concatenate((recoEmatch,recoPmatch)), ak.concatenate((genEmatch,genPmatch))

def getMatchesv2_lowPtOnly(tree,dRcut=0.3):
    ie = tree.GenEleClosestLpt.ind
    re = tree.GenEleClosestLpt.dr

    ip = tree.GenPosClosestLpt.ind
    rp = tree.GenPosClosestLpt.dr

    eMatch = (re<dRcut)
    ieMatch = ak.Array(ie[eMatch][:,np.newaxis].to_list())

    pMatch = (rp<dRcut)
    ipMatch = ak.Array(ip[pMatch][:,np.newaxis].to_list())

    electrons = tree.LptElectron
    
    genEmatch = tree.GenEle[eMatch]
    recoEmatch = electrons[eMatch][ieMatch][:,0]

    genPmatch = tree.GenPos[pMatch]
    recoPmatch = electrons[pMatch][ipMatch][:,0]
    
    return ak.concatenate((recoEmatch,recoPmatch)), ak.concatenate((genEmatch,genPmatch))

def getRegMatchesv2(tree,dRcut=0.3):
    ie = tree.GenEleClosest.ind
    te = tree.GenEleClosest.typ
    re = tree.GenEleClosest.dr

    ip = tree.GenPosClosest.ind
    tp = tree.GenPosClosest.typ
    rp = tree.GenPosClosest.dr

    eMatch = (ie>0) & (te==1) & (re<dRcut)
    ieMatch = ak.Array(ie[eMatch][:,np.newaxis].to_list())

    pMatch = (ip>0) & (tp==1) & (rp<dRcut)
    ipMatch = ak.Array(ip[pMatch][:,np.newaxis].to_list())

    genEmatch = tree.GenEle[eMatch]
    recoEmatch = tree.Electron[eMatch][ieMatch][:,0]

    genPmatch = tree.GenPos[pMatch]
    recoPmatch = tree.Electron[pMatch][ipMatch][:,0]
    
    return ak.concatenate((recoEmatch,recoPmatch)), ak.concatenate((genEmatch,genPmatch))

def getMatches(tree):
    eMatchLow = tree.GenEleMatch.typ==2
    pMatchLow = tree.GenPosMatch.typ==2
    eMatchInd = tree.GenEleMatch.ind
    pMatchInd = tree.GenPosMatch.ind
    sameMatch = eMatchLow & pMatchLow & (eMatchInd==pMatchInd)

    matchedGenEles = tree.GenEle[eMatchLow & ~sameMatch]
    matchedGenPos = tree.GenPos[pMatchLow & ~sameMatch]
    matchedGen_disamb = ak.where(tree.GenEleMatch.dr[sameMatch] < tree.GenPosMatch.dr[sameMatch],tree.GenEle[sameMatch],tree.GenPos[sameMatch])
    matchedInd_disamb = ak.Array(ak.where(tree.GenEleMatch.dr[sameMatch] < tree.GenPosMatch.dr[sameMatch],tree.GenEleMatch.ind[sameMatch],tree.GenPosMatch.ind[sameMatch])[:,np.newaxis].to_list())

    matchedToGenEle = tree.LptElectron[eMatchLow & ~sameMatch][ak.Array(tree.GenEleMatch[eMatchLow & ~sameMatch].ind[:,np.newaxis].to_list())][:,0]
    matchedToGenPos = tree.LptElectron[pMatchLow & ~sameMatch][ak.Array(tree.GenPosMatch[pMatchLow & ~sameMatch].ind[:,np.newaxis].to_list())][:,0]
    matchedToDisamb = tree.LptElectron[sameMatch][matchedInd_disamb][:,0]

    genMatched = ak.concatenate((matchedGenEles,matchedGenPos,matchedGen_disamb))
    
    LptMatchedPt = ak.concatenate((matchedToGenEle.pt,matchedToGenPos.pt,matchedToDisamb.pt))
    LptMatchedID = ak.concatenate((matchedToGenEle.ID,matchedToGenPos.ID,matchedToDisamb.ID))
    LptMatchedTrkProb = ak.concatenate((matchedToGenEle.trkProb,matchedToGenPos.trkProb,matchedToDisamb.trkProb))
    LptMatchedNumTrkHits = ak.concatenate((matchedToGenEle.numTrackerHits,matchedToGenPos.numTrackerHits,matchedToDisamb.numTrackerHits))
    
    # event numbers
    enE = tree.eventNum[eMatchLow & ~sameMatch]
    enP = tree.eventNum[pMatchLow & ~sameMatch]
    enD = tree.eventNum[sameMatch]
    en = ak.concatenate((enE,enP,enD))
    
    order = ak.argsort(en)
    
    return genMatched[order], LptMatchedPt[order], LptMatchedID[order], LptMatchedTrkProb[order], LptMatchedNumTrkHits[order]

def getRegMatches(tree):
    eMatch = tree.GenEleMatch.typ==1
    pMatch = tree.GenPosMatch.typ==1
    eMatchInd = tree.GenEleMatch.ind
    pMatchInd = tree.GenPosMatch.ind
    sameMatch = eMatch & pMatch & (eMatchInd==pMatchInd)

    matchedGenEles = tree.GenEle[eMatch & ~sameMatch]
    matchedGenPos = tree.GenPos[pMatch & ~sameMatch]
    matchedGen_disamb = ak.where(tree.GenEleMatch.dr[sameMatch] < tree.GenPosMatch.dr[sameMatch],tree.GenEle[sameMatch],tree.GenPos[sameMatch])
    matchedInd_disamb = ak.Array(ak.where(tree.GenEleMatch.dr[sameMatch] < tree.GenPosMatch.dr[sameMatch],tree.GenEleMatch.ind[sameMatch],tree.GenPosMatch.ind[sameMatch])[:,np.newaxis].to_list())

    matchedToGenEle = tree.Electron[eMatch & ~sameMatch][ak.Array(tree.GenEleMatch[eMatch & ~sameMatch].ind[:,np.newaxis].to_list())][:,0]
    matchedToGenPos = tree.Electron[pMatch & ~sameMatch][ak.Array(tree.GenPosMatch[pMatch & ~sameMatch].ind[:,np.newaxis].to_list())][:,0]
    matchedToDisamb = tree.Electron[sameMatch][matchedInd_disamb][:,0]

    genMatched = ak.concatenate((matchedGenEles,matchedGenPos,matchedGen_disamb))
    
    LptMatchedPt = ak.concatenate((matchedToGenEle.pt,matchedToGenPos.pt,matchedToDisamb.pt))
    #LptMatchedID = ak.concatenate((matchedToGenEle.ID,matchedToGenPos.ID,matchedToDisamb.ID))
    LptMatchedTrkProb = ak.concatenate((matchedToGenEle.trkProb,matchedToGenPos.trkProb,matchedToDisamb.trkProb))
    LptMatchedNumTrkHits = ak.concatenate((matchedToGenEle.numTrackerHits,matchedToGenPos.numTrackerHits,matchedToDisamb.numTrackerHits))
    
    # event numbers
    enE = tree.eventNum[eMatch & ~sameMatch]
    enP = tree.eventNum[pMatch & ~sameMatch]
    enD = tree.eventNum[sameMatch]
    en = ak.concatenate((enE,enP,enD))
    
    order = ak.argsort(en)
    
    return genMatched[order], LptMatchedPt[order], LptMatchedTrkProb[order], LptMatchedNumTrkHits[order]
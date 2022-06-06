#ifndef NTUPLECONTAINER_HH
#define NTUPLECONTAINER_HH

#include <vector>
using std::vector;
#include <iostream>

#include <TTree.h>

class NtupleContainer {

public:
    NtupleContainer();
    virtual ~NtupleContainer();
    void SetTree(TTree *tree);
    void CreateTreeBranches();
    void ClearTreeBranches();

    // Trigger and event-level branches
    unsigned long long eventNum_;
    unsigned long long runNum_;
    unsigned long long lumiSec_;
    bool isData_;

    // Gen branches
    
    // Gen particles
    int nGen_;
    float genwgt_;
    vector<int> genID_;
    // Only save hard-process gen particles
    vector<int> genCharge_;
    vector<float> genPt_;
    vector<float> genEta_;
    vector<float> genPhi_;
    vector<float> genEn_;
    vector<float> genPx_;
    vector<float> genPy_;
    vector<float> genPz_;
    vector<float> genVxy_;
    vector<float> genVz_;
    vector<float> genVx_;
    vector<float> genVy_;
    vector<float> genMass_;

    // Gen Electron & Positron from iDM signal
    int genEleCharge_;
    float genElePt_;
    float genEleEta_;
    float genElePhi_;
    float genEleEn_;
    float genElePx_;
    float genElePy_;
    float genElePz_;
    float genEleVxy_;
    float genEleVz_;
    float genEleVx_;
    float genEleVy_;
    int genEleBestMatchType_;
    int genEleBestMatchInd_;
    float genEleBestMatchDr_;
    int genEleBestMatchType_withCands_;
    int genEleBestMatchInd_withCands_;
    float genEleBestMatchDr_withCands_;
    int nGenEleMatches_;
    vector<int> genEleMatchTypes_;
    vector<int> genEleMatchInds_;
    vector<float> genEleMatchDrs_;

    int genPosCharge_;
    float genPosPt_;
    float genPosEta_;
    float genPosPhi_;
    float genPosEn_;
    float genPosPx_;
    float genPosPy_;
    float genPosPz_;
    float genPosVxy_;
    float genPosVz_;
    float genPosVx_;
    float genPosVy_;
    int genPosBestMatchType_;
    int genPosBestMatchInd_;
    float genPosBestMatchDr_;
    int genPosBestMatchType_withCands_;
    int genPosBestMatchInd_withCands_;
    float genPosBestMatchDr_withCands_;
    int nGenPosMatches_;
    vector<int> genPosMatchTypes_;
    vector<int> genPosMatchInds_;
    vector<float> genPosMatchDrs_;

    // Gen Electron + Positron info
    float genEEPt_;
    float genEEEta_;
    float genEEPhi_;
    float genEEEn_;
    float genEEMass_;
    float genEEdR_;
    
   // Gen jet
    int nGenJet_;
    vector<float> genJetPt_;
    vector<float> genJetEta_;
    vector<float> genJetPhi_;
    
    // Gen MET
    float genLeadMETPt_;
    float genLeadMETPhi_;
    float genLeadMETPx_;
    float genLeadMETPy_;
    float genLeadMETET_;

    // Reco Particles

    // Normal Electrons
    int nElectronDefault_;
    vector<float> recoElectronPt_;
    vector<float> recoElectronEta_;
    vector<float> recoElectronEtaError_;
    vector<float> recoElectronPhi_;
    vector<float> recoElectronPhiError_;
    vector<float> recoElectronAngularRes_;
    vector<float> recoElectronE_;
    vector<float> recoElectronPx_;
    vector<float> recoElectronPy_;
    vector<float> recoElectronPz_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<float> recoElectronDxy_;
    vector<float> recoElectronDxyError_;
    vector<float> recoElectronDz_;
    vector<float> recoElectronDzError_;
    vector<float> recoElectronTrkChi2_;
    vector<float> recoElectronTrkIso_;
    vector<int> recoElectronTrkNumTrackerHits_;
    vector<int> recoElectronTrkNumPixHits_;
    vector<int> recoElectronTrkNumStripHits_;
    vector<int> recoElectronCharge_;

    // Low pT electrons
    int nElectronLowPt_;
    vector<float> recoLowPtElectronPt_;
    vector<float> recoLowPtElectronEta_;
    vector<float> recoLowPtElectronEtaError_;
    vector<float> recoLowPtElectronPhi_;
    vector<float> recoLowPtElectronPhiError_;
    vector<float> recoLowPtElectronAngularRes_;
    vector<float> recoLowPtElectronE_;
    vector<float> recoLowPtElectronPx_;
    vector<float> recoLowPtElectronPy_;
    vector<float> recoLowPtElectronPz_;
    vector<float> recoLowPtElectronVxy_;
    vector<float> recoLowPtElectronVz_;
    vector<float> recoLowPtElectronDxy_;
    vector<float> recoLowPtElectronDxyError_;
    vector<float> recoLowPtElectronDz_;
    vector<float> recoLowPtElectronDzError_;
    vector<float> recoLowPtElectronTrkChi2_;
    vector<float> recoLowPtElectronTrkIso_;
    vector<int> recoLowPtElectronTrkNumTrackerHits_;
    vector<int> recoLowPtElectronTrkNumPixHits_;
    vector<int> recoLowPtElectronTrkNumStripHits_;
    vector<int> recoLowPtElectronCharge_;

    // Photons
    int nPhotons_;
    vector<float> PhotonEt_;
    vector<float> PhotonEta_;
    vector<float> PhotonPhi_;

    // OOT Photons
    int nOOTPhotons_;
    vector<float> ootPhotonEt_;
    vector<float> ootPhotonEta_;
    vector<float> ootPhotonPhi_;

    // Photon conversions
    int nConversions_;
    vector<float> conversionPt_;
    vector<float> conversionEta_;
    vector<float> conversionPhi_;
    vector<float> conversionE_;
    vector<float> conversionPx_;
    vector<float> conversionPy_;
    vector<float> conversionPz_;
    vector<float> conversionX_;
    vector<float> conversionY_;
    vector<float> conversionZ_;
    vector<float> conversionVxy_;
    vector<float> conversionVz_;
    vector<float> conversion_Trk1_innerPt_;
    vector<float> conversion_Trk1_innerEta_;
    vector<float> conversion_Trk1_innerPhi_;
    vector<float> conversion_Trk1_outerPt_;
    vector<float> conversion_Trk1_outerEta_;
    vector<float> conversion_Trk1_outerPhi_;
    vector<float> conversion_Trk2_innerPt_;
    vector<float> conversion_Trk2_innerEta_;
    vector<float> conversion_Trk2_innerPhi_;
    vector<float> conversion_Trk2_outerPt_;
    vector<float> conversion_Trk2_outerEta_;
    vector<float> conversion_Trk2_outerPhi_;

    // MET
    float PFMET_ET_;
    float PFMET_Px_;
    float PFMET_Py_;
    float PFMET_Pt_;
    float PFMET_Phi_;
    
    float CaloMET_ET_;
    float CaloMET_Px_;
    float CaloMET_Py_;
    float CaloMET_Pt_;
    float CaloMET_Phi_;

    float PuppiPFMET_ET_;
    float PuppiPFMET_Px_;
    float PuppiPFMET_Py_;
    float PuppiPFMET_Pt_;
    float PuppiPFMET_Phi_;
    
    float PuppiCaloMET_ET_;
    float PuppiCaloMET_Px_;
    float PuppiCaloMET_Py_;
    float PuppiCaloMET_Pt_;
    float PuppiCaloMET_Phi_;

    // Electron-positron vertices
    int nEleVertex_RR_;
    vector<int> RRvtx_idx1_;
    vector<int> RRvtx_idx2_;
    vector<float> RRvtx_recoVtxVxy_;
    vector<float> RRvtx_recoVtxSigmaVxy_;
    vector<float> RRvtx_recoVtxVx_;
    vector<float> RRvtx_recoVtxVy_;
    vector<float> RRvtx_recoVtxVz_;
    vector<float> RRvtx_recoVtxReducedChi2_;
    vector<float> RRvtx_prob_;
    vector<float> RRvtx_recoVtxDr_;
    vector<int> RRvtx_recoVtxSign_;
    vector<float> RRvtx_ll_pt_;
    vector<float> RRvtx_ll_eta_;
    vector<float> RRvtx_ll_phi_;
    vector<float> RRvtx_ll_e_;
    vector<float> RRvtx_ll_m_;
    vector<float> RRvtx_ll_px_;
    vector<float> RRvtx_ll_py_;
    vector<float> RRvtx_ll_pz_;

    int nEleVertex_LL_;
    vector<int> LLvtx_idx1_;
    vector<int> LLvtx_idx2_;
    vector<float> LLvtx_recoVtxVxy_;
    vector<float> LLvtx_recoVtxSigmaVxy_;
    vector<float> LLvtx_recoVtxVx_;
    vector<float> LLvtx_recoVtxVy_;
    vector<float> LLvtx_recoVtxVz_;
    vector<float> LLvtx_recoVtxReducedChi2_;
    vector<float> LLvtx_prob_;
    vector<float> LLvtx_recoVtxDr_;
    vector<int> LLvtx_recoVtxSign_;
    vector<float> LLvtx_ll_pt_;
    vector<float> LLvtx_ll_eta_;
    vector<float> LLvtx_ll_phi_;
    vector<float> LLvtx_ll_e_;
    vector<float> LLvtx_ll_m_;
    vector<float> LLvtx_ll_px_;
    vector<float> LLvtx_ll_py_;
    vector<float> LLvtx_ll_pz_;

    int nEleVertex_LR_;
    vector<int> LRvtx_idx1_;
    vector<int> LRvtx_idx2_;
    vector<float> LRvtx_recoVtxVxy_;
    vector<float> LRvtx_recoVtxSigmaVxy_;
    vector<float> LRvtx_recoVtxVx_;
    vector<float> LRvtx_recoVtxVy_;
    vector<float> LRvtx_recoVtxVz_;
    vector<float> LRvtx_recoVtxReducedChi2_;
    vector<float> LRvtx_prob_;
    vector<float> LRvtx_recoVtxDr_;
    vector<int> LRvtx_recoVtxSign_;
    vector<float> LRvtx_ll_pt_;
    vector<float> LRvtx_ll_eta_;
    vector<float> LRvtx_ll_phi_;
    vector<float> LRvtx_ll_e_;
    vector<float> LRvtx_ll_m_;
    vector<float> LRvtx_ll_px_;
    vector<float> LRvtx_ll_py_;
    vector<float> LRvtx_ll_pz_;

    int nEleVertex_RC_;
    vector<int> RCvtx_idx1_;
    vector<int> RCvtx_idx2_;
    vector<float> RCvtx_recoVtxVxy_;
    vector<float> RCvtx_recoVtxSigmaVxy_;
    vector<float> RCvtx_recoVtxVx_;
    vector<float> RCvtx_recoVtxVy_;
    vector<float> RCvtx_recoVtxVz_;
    vector<float> RCvtx_recoVtxReducedChi2_;
    vector<float> RCvtx_prob_;
    vector<float> RCvtx_recoVtxDr_;
    vector<int> RCvtx_recoVtxSign_;
    vector<float> RCvtx_ll_pt_;
    vector<float> RCvtx_ll_eta_;
    vector<float> RCvtx_ll_phi_;
    vector<float> RCvtx_ll_e_;
    vector<float> RCvtx_ll_m_;
    vector<float> RCvtx_ll_px_;
    vector<float> RCvtx_ll_py_;
    vector<float> RCvtx_ll_pz_;

    int nEleVertex_LC_;
    vector<int> LCvtx_idx1_;
    vector<int> LCvtx_idx2_;
    vector<float> LCvtx_recoVtxVxy_;
    vector<float> LCvtx_recoVtxSigmaVxy_;
    vector<float> LCvtx_recoVtxVx_;
    vector<float> LCvtx_recoVtxVy_;
    vector<float> LCvtx_recoVtxVz_;
    vector<float> LCvtx_recoVtxReducedChi2_;
    vector<float> LCvtx_prob_;
    vector<float> LCvtx_recoVtxDr_;
    vector<int> LCvtx_recoVtxSign_;
    vector<float> LCvtx_ll_pt_;
    vector<float> LCvtx_ll_eta_;
    vector<float> LCvtx_ll_phi_;
    vector<float> LCvtx_ll_e_;
    vector<float> LCvtx_ll_m_;
    vector<float> LCvtx_ll_px_;
    vector<float> LCvtx_ll_py_;
    vector<float> LCvtx_ll_pz_;

    // Electron candidates reconstructed from isoTracks + ECAL clusters
    int nEleCand_;
    vector<float> EleCand_pt_;
    vector<float> EleCand_et_;
    vector<float> EleCand_eta_;
    vector<float> EleCand_phi_;
    vector<float> EleCand_dxy_;
    vector<float> EleCand_dxyError_;
    vector<float> EleCand_dxy_PV_;
    vector<float> EleCand_dxyError_PV_;
    vector<float> EleCand_relPFiso_;
    vector<float> EleCand_relTrkiso_;
    vector<float> EleCand_ptDiff_;
    vector<float> EleCand_trkIso_;
    vector<float> EleCand_trkChi2_;
    vector<int> EleCand_numTrackerHits_;
    vector<int> EleCand_numPixHits_;
    vector<int> EleCand_numStripHits_;

    // Displaced dileptons from dedicated algorithm
    int ndispEE_;
    vector<int> dispEE_maxIxy_;
    vector<float> dispEE_Lxy_;
    vector<float> dispEE_Ixy_;
    vector<float> dispEE_trackDxy_;
    vector<float> dispEE_trackIxy_;
    vector<float> dispEE_vx_;
    vector<float> dispEE_vy_;
    vector<float> dispEE_mass_;
    vector<float> dispEE_normalizedChi2_;
    vector<float> dispEE_leadingPt_;
    vector<float> dispEE_subleadingPt_;
    vector<float> dispEE_leadingEt_;
    vector<float> dispEE_subleadingEt_;
    vector<float> dispEE_cosAlpha_;
    vector<float> dispEE_dPhi_;
    vector<float> dispEE_relisoA_;
    vector<float> dispEE_relisoB_;
    vector<int> dispEE_fromPVA_;
    vector<int> dispEE_fromPVB_;
    vector<int> dispEE_PVAssociation_;

    // Displaced dilepton candidates (aren't required to pass a baseline selection)
    int nEECand_;
    vector<float> EECand_Lxy_PV_;
    vector<float> EECand_Ixy_PV_;
    vector<float> EECand_Lxy_0_;
    vector<float> EECand_Ixy_0_;
    vector<float> EECand_Lxy_BS_;
    vector<float> EECand_Ixy_BS_;
    vector<float> EECand_trackDxy_;
    vector<float> EECand_trackIxy_;
    vector<float> EECand_trackDxy_PV_;
    vector<float> EECand_trackIxy_PV_;
    vector<float> EECand_trackDxy_0_;
    vector<float> EECand_trackIxy_0_;
    vector<float> EECand_trackDxy_BS_;
    vector<float> EECand_trackIxy_BS_;
    vector<float> EECand_vx_;
    vector<float> EECand_vy_;
    vector<float> EECand_mass_;
    vector<float> EECand_normalizedChi2_;
    vector<float> EECand_leadingPt_;
    vector<float> EECand_subleadingPt_;
    vector<float> EECand_leadingEt_;
    vector<float> EECand_subleadingEt_;
    vector<float> EECand_cosAlpha_;
    vector<float> EECand_dR_;
    vector<float> EECand_dPhi_;
    vector<float> EECand_lldPhi_;
    vector<float> EECand_relisoA_;
    vector<float> EECand_relisoB_;

private:
    // Reco and gen TTrees
    TTree * outT;

};


#endif
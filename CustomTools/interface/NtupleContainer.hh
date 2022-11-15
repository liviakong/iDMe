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
    unsigned int fired_;
    unsigned int fired16_;
    unsigned int fired17_;
    unsigned int fired18_;
    unsigned long long eventNum_;
    unsigned long long runNum_;
    unsigned long long lumiSec_;
    bool isData_;

    // MET Filters
    uint32_t METFiltersFailBits_;

    //////////////////////
    //// Gen branches ////
    //////////////////////
    
    // Gen particles
    int nGen_;
    float genwgt_;
    int genpuobs_;
    int genputrue_;
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
    int genEleMotherID_;
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
    int genEleClosestType_;
    int genEleClosestInd_;
    float genEleClosestDr_;
    int genEleClosestInd_reg_;
    float genEleClosestDr_reg_;
    int genEleClosestInd_lpt_;
    float genEleClosestDr_lpt_;


    int genPosCharge_;
    int genPosMotherID_;
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
    int genPosClosestType_;
    int genPosClosestInd_;
    float genPosClosestDr_;
    int genPosClosestInd_reg_;
    float genPosClosestDr_reg_;
    int genPosClosestInd_lpt_;
    float genPosClosestDr_lpt_;

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
    vector<bool> recoElectronMatched_;
    vector<int> recoElectronMatchType_;
    vector<float> recoElectronMatchdR_;
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
    vector<float> recoElectronTrkRelIso_;
    vector<float> recoElectronTrkProb_;
    vector<int> recoElectronTrkNumTrackerHits_;
    vector<int> recoElectronTrkNumPixHits_;
    vector<int> recoElectronTrkNumStripHits_;
    vector<int> recoElectronCharge_;

    // Low pT electrons
    int nElectronLowPt_;
    vector<bool> recoLowPtElectronMatched_;
    vector<int> recoLowPtElectronMatchType_;
    vector<float> recoLowPtElectronMatchdR_;
    vector<float> recoLowPtElectronPt_;
    vector<float> recoLowPtElectronEta_;
    vector<float> recoLowPtElectronEtaError_;
    vector<float> recoLowPtElectronPhi_;
    vector<float> recoLowPtElectronPhiError_;
    vector<float> recoLowPtElectronID_;
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
    vector<float> recoLowPtElectronTrkRelIso_;
    vector<float> recoLowPtElectronTrkProb_;
    vector<int> recoLowPtElectronTrkNumTrackerHits_;
    vector<int> recoLowPtElectronTrkNumPixHits_;
    vector<int> recoLowPtElectronTrkNumStripHits_;
    vector<int> recoLowPtElectronCharge_;
    vector<float> recoLowPtElectronMinDrToReg_;

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

    // Jets
    int PFNJet_;
    int PFNPassIDJet_;
    int PFNHighPtJet_;
    vector<float> PFJetPt_;
    vector<float> PFJetEta_;
    vector<float> PFJetPhi_;
    vector<float> PFJetCorrectedPt_;
    vector<float> PFJetCorrectedEta_;
    vector<float> PFJetCorrectedPhi_;
    vector<float> PFJetCorrectedBTag_;
    vector<float> PFJetCorrectedCHEF_;
    vector<float> PFJetCorrectedNHEF_;
    vector<float> PFJetCorrectedCEEF_;
    vector<float> PFJetCorrectedNEEF_;
    vector<float> PFJetCorrectedNumDaughters_;
    vector<float> PFJetCorrectedChargedMultiplicity_;
    vector<float> PFJetCorrectedJESUpPt_;
    vector<float> PFJetCorrectedJESUpEta_;
    vector<float> PFJetCorrectedJESUpPhi_;
    vector<float> PFJetCorrectedJESDownPt_;
    vector<float> PFJetCorrectedJESDownEta_;
    vector<float> PFJetCorrectedJESDownPhi_;
    vector<float> PFJetCorrectedJERUpPt_;
    vector<float> PFJetCorrectedJERUpEta_;
    vector<float> PFJetCorrectedJERUpPhi_;
    vector<float> PFJetCorrectedJERDownPt_;
    vector<float> PFJetCorrectedJERDownEta_;
    vector<float> PFJetCorrectedJERDownPhi_;
    bool PFHEMFlag_;

    // MET
    float PFMET_ET_;
    float PFMET_Px_;
    float PFMET_Py_;
    float PFMET_Pt_;
    float PFMET_Phi_;
    float PFMETSmearingOnlyPt_;
    float PFMETSmearingOnlyPhi_;
    float PFMETCorrectedPt_;
    float PFMETCorrectedPhi_;
    float PFMETEEDeltaPx_;
    float PFMETEEDeltaPy_;
    float PFMETJESUpPt_;
    float PFMETJESUpPhi_;
    float PFMETJESDownPt_;
    float PFMETJESDownPhi_;
    float PFMETJERUpPt_;
    float PFMETJERUpPhi_;
    float PFMETJERDownPt_;
    float PFMETJERDownPhi_;
    float PFMETMuonEtFraction_;
    
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

    // Pileup density
    float rho_;

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

protected:
    // Reco and gen TTrees
    TTree * outT;

};


#endif
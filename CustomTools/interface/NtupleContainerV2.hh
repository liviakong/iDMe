#ifndef NTUPLECONTAINERV2_HH
#define NTUPLECONTAINERV2_HH

#include <vector>
using std::vector;
#include <iostream>

#include <TTree.h>

class NtupleContainerV2 {

public:
    NtupleContainerV2();
    virtual ~NtupleContainerV2();
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
    bool isSignal_;

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
    int genPosClosestType_;
    int genPosClosestInd_;
    float genPosClosestDr_;
    int genPosClosestInd_reg_;
    float genPosClosestDr_reg_;
    int genPosClosestInd_lpt_;
    float genPosClosestDr_lpt_;

    // Gen Electron matches to isoTracks
    int nGenEleTrkMatches;
    float genEleNearestTrack_pt;
    float genEleNearestTrack_eta;
    float genEleNearestTrack_phi;
    float genEleNearestTrack_dRGen;
    float genEleNearestTrack_pfIso3;
    float genEleNearestTrack_miniIso;
    float genEleNearestTrack_dxy;
    float genEleNearestTrack_dz;
    bool genEleNearestTrack_highPurity;
    bool genEleNearestTrack_Loose;
    int genEleNearestTrack_charge;
    int genEleNearestTrack_numPixHits;
    int genEleNearestTrack_numStripHits;
    int genEleNearestTrack_fromPV;
    int genEleNearestTrack_tkIdx;

    // Gen Positron matches to isoTracks
    int nGenPosTrkMatches;
    float genPosNearestTrack_pt;
    float genPosNearestTrack_eta;
    float genPosNearestTrack_phi;
    float genPosNearestTrack_dRGen;
    float genPosNearestTrack_pfIso3;
    float genPosNearestTrack_miniIso;
    float genPosNearestTrack_dxy;
    float genPosNearestTrack_dz;
    bool genPosNearestTrack_highPurity;
    bool genPosNearestTrack_Loose;
    int genPosNearestTrack_charge;
    int genPosNearestTrack_numPixHits;
    int genPosNearestTrack_numStripHits;
    int genPosNearestTrack_fromPV;
    int genPosNearestTrack_tkIdx;

    // Gen Electron + Positron info
    float genEEPt_;
    float genEEEta_;
    float genEEPhi_;
    float genEEEn_;
    float genEEMass_;
    float genEEdR_;
    float genEEMETdPhi_;
    
   // Gen jet
    int nGenJet_;
    vector<float> genJetPt_;
    vector<float> genJetEta_;
    vector<float> genJetPhi_;
    vector<float> genJetMETdPhi_;
    
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
    vector<float> recoElectronID_cutVeto_;
    vector<float> recoElectronID_cutLoose_;
    vector<float> recoElectronID_cutMed_;
    vector<float> recoElectronID_cutTight_;
    vector<float> recoElectronID_mvaIso90_;
    vector<float> recoElectronID_mvaIso80_;
    vector<float> recoElectronID_mvaIsoLoose_;
    vector<float> recoElectronID_mva90_;
    vector<float> recoElectronID_mva80_;
    vector<float> recoElectronID_mvaLoose_;
    vector<float> recoElectronAngularRes_;
    vector<float> recoElectronE_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<float> recoElectronDxy_;
    vector<float> recoElectronDxyError_;
    vector<float> recoElectronDz_;
    vector<float> recoElectronDzError_;
    vector<float> recoElectronTrkChi2_;
    vector<float> recoElectronTrkIso_;
    vector<float> recoElectronTrkRelIso_;
    vector<float> recoElectronCaloIso_;
    vector<float> recoElectronCaloRelIso_;
    vector<float> recoElectronPFIso_dR4_;
    vector<float> recoElectronPFRelIso_dR4_;
    vector<float> recoElectronPFIso_dR3_;
    vector<float> recoElectronPFRelIso_dR3_;
    vector<float> recoElectronPFIso_dR8_;
    vector<float> recoElectronPFRelIso_dR8_;
    vector<float> recoElectronPFIso_;
    vector<float> recoElectronPFRelIso_;
    vector<float> recoElectronChadIso_;
    vector<float> recoElectronNhadIso_;
    vector<float> recoElectronPhoIso_;
    vector<float> recoElectronRhoEA_;
    vector<float> recoElectronTrkProb_;
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
    vector<float> recoLowPtElectronID_;
    vector<float> recoLowPtElectronAngularRes_;
    vector<float> recoLowPtElectronE_;
    vector<float> recoLowPtElectronVxy_;
    vector<float> recoLowPtElectronVz_;
    vector<float> recoLowPtElectronDxy_;
    vector<float> recoLowPtElectronDxyError_;
    vector<float> recoLowPtElectronDz_;
    vector<float> recoLowPtElectronDzError_;
    vector<float> recoLowPtElectronTrkChi2_;
    vector<float> recoLowPtElectronTrkIso_;
    vector<float> recoLowPtElectronTrkRelIso_;
    vector<float> recoLowPtElectronCaloIso_;
    vector<float> recoLowPtElectronCaloRelIso_;
    vector<float> recoLowPtElectronPFIso_dR4_;
    vector<float> recoLowPtElectronPFRelIso_dR4_;
    vector<float> recoLowPtElectronPFIso_dR3_;
    vector<float> recoLowPtElectronPFRelIso_dR3_;
    vector<float> recoLowPtElectronPFIso_dR8_;
    vector<float> recoLowPtElectronPFRelIso_dR8_;
    vector<float> recoLowPtElectronTrkProb_;
    vector<float> recoLowPtElectronPFIso_;
    vector<float> recoLowPtElectronPFRelIso_;
    vector<float> recoLowPtElectronChadIso_;
    vector<float> recoLowPtElectronNhadIso_;
    vector<float> recoLowPtElectronPhoIso_;
    vector<float> recoLowPtElectronRhoEA_;
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
    int PFNJetAll_;
    vector<float> PFJetPt_;
    vector<float> PFJetEta_;
    vector<float> PFJetPhi_;
    vector<float> PFJetBTag_;
    vector<float> PFJetMETdPhi_;
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
    float PFMET_Pt_;
    float PFMET_Phi_;
    float PFMETJESUpPt_;
    float PFMETJESUpPhi_;
    float PFMETJESDownPt_;
    float PFMETJESDownPhi_;
    float PFMETJERUpPt_;
    float PFMETJERUpPhi_;
    float PFMETJERDownPt_;
    float PFMETJERDownPhi_;
    float PFMETJetResUpSmearPt_;
    float PFMETJetResUpSmearPhi_;
    float PFMETJetResDownSmearPt_;
    float PFMETJetResDownSmearPhi_;
    
    float CaloMET_ET_;
    float CaloMET_Pt_;
    float CaloMET_Phi_;

    // Pileup density
    float rho_;

    // Electron-positron vertices
    int nvtx_;
    vector<std::string> vtx_type_;
    vector<float> vtx_recoVtxVxy_;
    vector<float> vtx_recoVtxSigmaVxy_;
    vector<float> vtx_recoVtxVx_;
    vector<float> vtx_recoVtxVy_;
    vector<float> vtx_recoVtxVz_;
    vector<float> vtx_recoVtxReducedChi2_;
    vector<float> vtx_prob_;
    vector<float> vtx_recoVtxDr_;
    vector<int> vtx_recoVtxSign_;
    vector<float> vtx_minDxy_;
    vector<float> vtx_METdPhi_;
    vector<float> vtx_ll_pt_;
    vector<float> vtx_ll_eta_;
    vector<float> vtx_ll_phi_;
    vector<float> vtx_ll_e_;
    vector<float> vtx_ll_m_;
    vector<float> vtx_ll_px_;
    vector<float> vtx_ll_py_;
    vector<float> vtx_ll_pz_;
    vector<float> vtx_ll_PFIso_dR4_;
    vector<float> vtx_ll_PFRelIso_dR4_;
    vector<float> vtx_ll_PFIso_dR3_;
    vector<float> vtx_ll_PFRelIso_dR3_;
    vector<float> vtx_ll_PFIso_dR8_;
    vector<float> vtx_ll_PFRelIso_dR8_;

    vector<std::string> vtx_e1_type_;
    vector<int> vtx_e1_idx_;
    vector<std::string> vtx_e2_type_;
    vector<int> vtx_e2_idx_;

protected:
    // Reco and gen TTrees
    TTree * outT;

};


#endif
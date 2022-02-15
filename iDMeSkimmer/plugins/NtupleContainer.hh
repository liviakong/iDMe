#ifndef NTUPLECONTAINER_HH
#define NTUPLECONTAINER_HH

#include <vector>
using std::vector;

#include <TTree.h>

class NtupleContainer {

public:
    NtupleContainer();
    virtual ~NtupleContainer();
    void SetTree(TTree *tree, bool isData);
    void CreateTreeBranches();
    void ClearTreeBranches();

    // Trigger and event-level branches
    unsigned long long eventNum_;
    unsigned long long runNum_;
    unsigned long long lumiSec_;

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
    vector<float> genVtx_x_;
    vector<float> genVtx_y_;
    vector<float> genVtx_z_;
    vector<float> genMass_;
    
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
    vector<float> recoElectronPhi_;
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
    vector<int> recoElectron_passConversionVeto_;

    // Low pT electrons
    int nElectronLowPt_;
    vector<float> recoLowPtElectronPt_;
    vector<float> recoLowPtElectronEta_;
    vector<float> recoLowPtElectronPhi_;
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
    vector<int> recoLowPtElectron_passConversionVeto_;

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
    int nEleVertex_regreg_;
    vector<float> regreg_recoVtxVxy_;
    vector<float> regreg_recoVtxVz_;
    vector<float> regreg_recoVtxSigmaVxy_;
    vector<float> regreg_recoVtx_x_;
    vector<float> regreg_recoVtx_y_;
    vector<float> regreg_recoVtx_z_;
    vector<float> regreg_recoVtxReducedChi2_;
    vector<float> regreg_recoVtxDr_;
    vector<int> regreg_recoVtxSign_;
    vector<float> regreg_ll_pt_;
    vector<float> regreg_ll_eta_;
    vector<float> regreg_ll_phi_;
    vector<float> regreg_ll_e_;
    vector<float> regreg_ll_m_;
    vector<float> regreg_ll_px_;
    vector<float> regreg_ll_py_;
    vector<float> regreg_ll_pz_;

    int nEleVertex_lowlow_;
    vector<float> lowlow_recoVtxVxy_;
    vector<float> lowlow_recoVtxVz_;
    vector<float> lowlow_recoVtxSigmaVxy_;
    vector<float> lowlow_recoVtx_x_;
    vector<float> lowlow_recoVtx_y_;
    vector<float> lowlow_recoVtx_z_;
    vector<float> lowlow_recoVtxReducedChi2_;
    vector<float> lowlow_recoVtxDr_;
    vector<int> lowlow_recoVtxSign_;
    vector<float> lowlow_ll_pt_;
    vector<float> lowlow_ll_eta_;
    vector<float> lowlow_ll_phi_;
    vector<float> lowlow_ll_e_;
    vector<float> lowlow_ll_m_;
    vector<float> lowlow_ll_px_;
    vector<float> lowlow_ll_py_;
    vector<float> lowlow_ll_pz_;

    int nEleVertex_lowreg_;
    vector<float> lowreg_recoVtxVxy_;
    vector<float> lowreg_recoVtxVz_;
    vector<float> lowreg_recoVtxSigmaVxy_;
    vector<float> lowreg_recoVtx_x_;
    vector<float> lowreg_recoVtx_y_;
    vector<float> lowreg_recoVtx_z_;
    vector<float> lowreg_recoVtxReducedChi2_;
    vector<float> lowreg_recoVtxDr_;
    vector<int> lowreg_recoVtxSign_;
    vector<float> lowreg_ll_pt_;
    vector<float> lowreg_ll_eta_;
    vector<float> lowreg_ll_phi_;
    vector<float> lowreg_ll_e_;
    vector<float> lowreg_ll_m_;
    vector<float> lowreg_ll_px_;
    vector<float> lowreg_ll_py_;
    vector<float> lowreg_ll_pz_;

private:
    // Reco and gen TTrees
    TTree * outT;
    bool isData_;

};


#endif
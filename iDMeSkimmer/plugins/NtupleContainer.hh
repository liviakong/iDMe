#ifndef NTUPLECONTAINER_HH
#define NTUPLECONTAINER_HH

#include <vector>
using std::vector;

#include <TTree.h>

class NtupleContainer {

public:
    NtupleContainer();
    virtual ~NtupleContainer();
    void SetRecoTree(TTree *tree);
    void SetGenTree(TTree *tree);
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
    //std::vector<bool> genHardProcess_; --> would always be 1 
    vector<int> genCharge_;
    vector<float> genPt_;
    vector<float> genEta_;
    vector<float> genPhi_;
    vector<float> genPz_;
    vector<float> genEn_;
    vector<float> genVxy_;
    vector<float> genVz_;
    vector<float> genMass_;
    
   // Gen jet
    vector<float> genJetPt_;
    vector<float> genJetEta_;
    vector<float> genJetPhi_;
    
    // Gen MET
    float genLeadMETPt_;
    float genLeadMETPhi_;

    // Reco Particles

    // Normal Electrons
    int nElectronDefault_;
    vector<float> recoElectronPt_;
    vector<float> recoElectronEta_;
    vector<float> recoElectronPhi_;
    vector<float> recoElectronVxy_;
    vector<float> recoElectronVz_;
    vector<int> recoElectronCharge_;

    // Low pT electrons
    int nElectronLowPt_;
    vector<float> recoLowPtElectronPt_;
    vector<float> recoLowPtElectronEta_;
    vector<float> recoLowPtElectronPhi_;
    vector<float> recoLowPtElectronVxy_;
    vector<float> recoLowPtElectronVz_;
    vector<int> recoLowPtElectronCharge_;

    // Electron-positron vertices
    vector<float> regreg_recoVtxVxy_;
    vector<float> regreg_recoVtxVz_;
    vector<float> regreg_recoVtxSigmaVxy_;
    vector<float> regreg_recoVtxReducedChi2_;
    vector<float> regreg_recoVtxDr_;
    vector<float> lowlow_recoVtxVxy_;
    vector<float> lowlow_recoVtxVz_;
    vector<float> lowlow_recoVtxSigmaVxy_;
    vector<float> lowlow_recoVtxReducedChi2_;
    vector<float> lowlow_recoVtxDr_;
    vector<float> lowreg_recoVtxVxy_;
    vector<float> lowreg_recoVtxVz_;
    vector<float> lowreg_recoVtxSigmaVxy_;
    vector<float> lowreg_recoVtxReducedChi2_;
    vector<float> lowreg_recoVtxDr_;

private:
    // Reco and gen TTrees
    TTree * recoT;
    TTree * genT;
    bool isData_;

};


#endif
#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetRecoTree(TTree *tree) { recoT = tree; }
void NtupleContainer::SetGenTree(TTree *tree) { genT = tree; isData_ = false; }

void NtupleContainer::CreateTreeBranches() {

    // Reco information
    recoT->Branch("eventNum", &eventNum_);
    recoT->Branch("lumiSec", &lumiSec_);
    recoT->Branch("runNum", &runNum_);

    // Normal Electrons
    recoT->Branch("nElectron",&nElectronDefault_);
    recoT->Branch("Electron_pt",&recoElectronPt_);
    recoT->Branch("Electron_eta",&recoElectronEta_);
    recoT->Branch("Electron_phi",&recoElectronPhi_);
    recoT->Branch("Electron_vxy",&recoElectronVxy_);
    recoT->Branch("Electron_vz",&recoElectronVz_);
    recoT->Branch("Electron_dxy",&recoElectronDxy_);
    recoT->Branch("Electron_dxyErr",&recoElectronDxyError_);
    recoT->Branch("Electron_dz",&recoElectronDz_);
    recoT->Branch("Electron_dzErr",&recoElectronDzError_);
    recoT->Branch("Electron_trkChi2",&recoElectronTrkChi2_);
    recoT->Branch("Electron_trkIso",&recoElectronTrkIso_);
    recoT->Branch("Electron_numTrackerHits",&recoElectronTrkNumTrackerHits_);
    recoT->Branch("Electron_numPixHits",&recoElectronTrkNumPixHits_);
    recoT->Branch("Electron_numStripHits",&recoElectronTrkNumStripHits_);
    recoT->Branch("Electron_charge",&recoElectronCharge_);

    // Low pT electrons
    recoT->Branch("nLptElectron",&nElectronLowPt_);
    recoT->Branch("LptElectron_pt",&recoLowPtElectronPt_);
    recoT->Branch("LptElectron_eta",&recoLowPtElectronEta_);
    recoT->Branch("LptElectron_phi",&recoLowPtElectronPhi_);
    recoT->Branch("LptElectron_vxy",&recoLowPtElectronVxy_);
    recoT->Branch("LptElectron_vz",&recoLowPtElectronVz_);
    recoT->Branch("LptElectron_dxy",&recoLowPtElectronDxy_);
    recoT->Branch("LptElectron_dxyErr",&recoLowPtElectronDxyError_);
    recoT->Branch("LptElectron_dz",&recoLowPtElectronDz_);
    recoT->Branch("LptElectron_dzErr",&recoLowPtElectronDzError_);
    recoT->Branch("LptElectron_trkChi2",&recoLowPtElectronTrkChi2_);
    recoT->Branch("LptElectron_trkIso",&recoLowPtElectronTrkIso_);
    recoT->Branch("LptElectron_numTrackerHits",&recoLowPtElectronTrkNumTrackerHits_);
    recoT->Branch("LptElectron_numPixHits",&recoLowPtElectronTrkNumPixHits_);
    recoT->Branch("LptElectron_numStripHits",&recoLowPtElectronTrkNumStripHits_);
    recoT->Branch("LptElectron_charge",&recoLowPtElectronCharge_);

    // Electron-positron vertex branches
    recoT->Branch("nEleVertex_regreg",&nEleVertex_regreg_);
    recoT->Branch("EleVertex_regreg_vxy", &regreg_recoVtxVxy_);
    recoT->Branch("EleVertex_regreg_vz",  &regreg_recoVtxVz_);
    recoT->Branch("EleVertex_regreg_sigmavxy", &regreg_recoVtxSigmaVxy_);
    recoT->Branch("EleVertex_regreg_reduced_chi2", &regreg_recoVtxReducedChi2_);
    recoT->Branch("EleVertex_regreg_dR",  &regreg_recoVtxDr_);
    recoT->Branch("nEleVertex_lowlow",&nEleVertex_lowlow_);
    recoT->Branch("EleVertex_lowlow_vxy", &lowlow_recoVtxVxy_);
    recoT->Branch("EleVertex_lowlow_vz",  &lowlow_recoVtxVz_);
    recoT->Branch("EleVertex_lowlow_sigmavxy", &lowlow_recoVtxSigmaVxy_);
    recoT->Branch("EleVertex_lowlow_reduced_chi2", &lowlow_recoVtxReducedChi2_);
    recoT->Branch("EleVertex_lowlow_dR",  &lowlow_recoVtxDr_);
    recoT->Branch("nEleVertex_lowreg",&nEleVertex_lowreg_);
    recoT->Branch("EleVertex_lowreg_vxy", &lowreg_recoVtxVxy_);
    recoT->Branch("EleVertex_lowreg_vz",  &lowreg_recoVtxVz_);
    recoT->Branch("EleVertex_lowreg_sigmavxy", &lowreg_recoVtxSigmaVxy_);
    recoT->Branch("EleVertex_lowreg_reduced_chi2", &lowreg_recoVtxReducedChi2_);
    recoT->Branch("EleVertex_lowreg_dR",  &lowreg_recoVtxDr_);

    // Gen information
    if (!isData_) {
        genT->Branch("eventNum", &eventNum_);
        genT->Branch("genWgt", &genwgt_);
        genT->Branch("nGenPart",&nGen_);
        genT->Branch("GenPart_ID", &genID_);
        genT->Branch("GenPart_charge", &genCharge_);
        genT->Branch("GenPart_pt", &genPt_);
        genT->Branch("GenPart_eta", &genEta_);
        genT->Branch("GenPart_phi", &genPhi_);
        genT->Branch("GenPart_pz", &genPz_);
        genT->Branch("GenPart_energy", &genEn_);
        genT->Branch("GenPart_vxy", &genVxy_);
        genT->Branch("GenPart_vz", &genVz_);
        genT->Branch("GenPart_mass", &genMass_);
        genT->Branch("nGenJet",&nGenJet_);
        genT->Branch("GenJet_pt", &genJetPt_);
        genT->Branch("GenJet_eta", &genJetEta_);
        genT->Branch("GenJet_phi", &genJetPhi_);
        genT->Branch("GenMET_pt", &genLeadMETPt_);
        genT->Branch("GenMET_phi", &genLeadMETPhi_);
    }

}

void NtupleContainer::ClearTreeBranches() {
    
    // Gen particles
    nGen_ = 0;
    genID_.clear();
    genCharge_.clear();
    genPt_.clear();
    genEta_.clear();
    genPhi_.clear();
    genPz_.clear();
    genEn_.clear();
    genVxy_.clear();
    genVz_.clear();
    genMass_.clear();
    
    // Gen jet
    nGenJet_ = 0;
    genJetPt_.clear();
    genJetEta_.clear();
    genJetPhi_.clear();

    //Gen MET
    genLeadMETPt_ = -9999;
    genLeadMETPhi_ = -9999;

    // Electrons
    nElectronDefault_ = 0;
    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronDxy_.clear();
    recoElectronDxyError_.clear();
    recoElectronDz_.clear();
    recoElectronDzError_.clear();
    recoElectronTrkChi2_.clear();
    recoElectronTrkIso_.clear();
    recoElectronTrkNumTrackerHits_.clear();
    recoElectronTrkNumPixHits_.clear();
    recoElectronTrkNumStripHits_.clear();
    recoElectronCharge_.clear();

    // Low pT electrons
    nElectronLowPt_ = 0;
    recoLowPtElectronPt_.clear();
    recoLowPtElectronEta_.clear();
    recoLowPtElectronPhi_.clear();
    recoLowPtElectronVxy_.clear();
    recoLowPtElectronVz_.clear();
    recoLowPtElectronDxy_.clear();
    recoLowPtElectronDxyError_.clear();
    recoLowPtElectronDz_.clear();
    recoLowPtElectronDzError_.clear();
    recoLowPtElectronTrkChi2_.clear();
    recoLowPtElectronTrkIso_.clear();
    recoLowPtElectronTrkNumTrackerHits_.clear();
    recoLowPtElectronTrkNumPixHits_.clear();
    recoLowPtElectronTrkNumStripHits_.clear();
    recoLowPtElectronCharge_.clear();

    // Electron-positron vertices
    nEleVertex_regreg_ = 0;
    regreg_recoVtxVxy_.clear();
    regreg_recoVtxVz_.clear();
    regreg_recoVtxSigmaVxy_.clear();
    regreg_recoVtxReducedChi2_.clear();
    regreg_recoVtxDr_.clear();
    nEleVertex_lowlow_ = 0;
    lowlow_recoVtxVxy_.clear();
    lowlow_recoVtxVz_.clear();
    lowlow_recoVtxSigmaVxy_.clear();
    lowlow_recoVtxReducedChi2_.clear();
    lowlow_recoVtxDr_.clear();
    nEleVertex_lowreg_ = 0;
    lowreg_recoVtxVxy_.clear();
    lowreg_recoVtxVz_.clear();
    lowreg_recoVtxSigmaVxy_.clear();
    lowreg_recoVtxReducedChi2_.clear();
    lowreg_recoVtxDr_.clear();
}
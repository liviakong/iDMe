#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetRecoTree(TTree *tree) { recoT = tree; }
void NtupleContainer::SetGenTree(TTree *tree) { genT = tree; isData_ = false; }

void NtupleContainer::CreateTreeBranches() {

    // Reco information
    recoT->Branch("event_num", &eventNum_);
    recoT->Branch("lumi_sec", &lumiSec_);
    recoT->Branch("run_num", &runNum_);

    // Normal Electrons
    recoT->Branch("n_ele",&nElectronDefault_);
    recoT->Branch("ele_pt",&recoElectronPt_);
    recoT->Branch("ele_eta",&recoElectronEta_);
    recoT->Branch("ele_phi",&recoElectronPhi_);
    recoT->Branch("ele_vxy",&recoElectronVxy_);
    recoT->Branch("ele_vz",&recoElectronVz_);
    recoT->Branch("ele_charge",&recoElectronCharge_);

    // Low pT electrons
    recoT->Branch("n_lowpt_ele",&nElectronLowPt_);
    recoT->Branch("lowpt_ele_pt",&recoLowPtElectronPt_);
    recoT->Branch("lowpt_ele_eta",&recoLowPtElectronEta_);
    recoT->Branch("lowpt_ele_phi",&recoLowPtElectronPhi_);
    recoT->Branch("lowpt_ele_vxy",&recoLowPtElectronVxy_);
    recoT->Branch("lowpt_ele_vz",&recoLowPtElectronVz_);
    recoT->Branch("lowpt_ele_charge",&recoLowPtElectronCharge_);

    // Electron-positron vertex branches
    recoT->Branch("reco_vtx_regreg_vxy", &regreg_recoVtxVxy_);
    recoT->Branch("reco_vtx_regreg_vz",  &regreg_recoVtxVz_);
    recoT->Branch("reco_vtx_regreg_sigmavxy", &regreg_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_regreg_reduced_chi2", &regreg_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_regreg_dR",  &regreg_recoVtxDr_);
    recoT->Branch("reco_vtx_lowlow_vxy", &lowlow_recoVtxVxy_);
    recoT->Branch("reco_vtx_lowlow_vz",  &lowlow_recoVtxVz_);
    recoT->Branch("reco_vtx_lowlow_sigmavxy", &lowlow_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_lowlow_reduced_chi2", &lowlow_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_lowlow_dR",  &lowlow_recoVtxDr_);
    recoT->Branch("reco_vtx_lowreg_vxy", &lowreg_recoVtxVxy_);
    recoT->Branch("reco_vtx_lowreg_vz",  &lowreg_recoVtxVz_);
    recoT->Branch("reco_vtx_lowreg_sigmavxy", &lowreg_recoVtxSigmaVxy_);
    recoT->Branch("reco_vtx_lowreg_reduced_chi2", &lowreg_recoVtxReducedChi2_);
    recoT->Branch("reco_vtx_lowreg_dR",  &lowreg_recoVtxDr_);

    // Gen information
    if (!isData_) {
        genT->Branch("event_num", &eventNum_);
        genT->Branch("gen_wgt", &genwgt_);
        genT->Branch("n_gen",&nGen_);
        genT->Branch("gen_ID", &genID_);
        genT->Branch("gen_charge", &genCharge_);
        genT->Branch("gen_pt", &genPt_);
        genT->Branch("gen_eta", &genEta_);
        genT->Branch("gen_phi", &genPhi_);
        genT->Branch("gen_pz", &genPz_);
        genT->Branch("gen_energy", &genEn_);
        genT->Branch("gen_vxy", &genVxy_);
        genT->Branch("gen_vz", &genVz_);
        genT->Branch("gen_mass", &genMass_);
        genT->Branch("gen_jet_pt", &genJetPt_);
        genT->Branch("gen_jet_eta", &genJetEta_);
        genT->Branch("gen_jet_phi", &genJetPhi_);
        genT->Branch("gen_MET_pt", &genLeadMETPt_);
        genT->Branch("gen_MET_phi", &genLeadMETPhi_);
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
    recoElectronCharge_.clear();

    // Low pT electrons
    nElectronLowPt_ = 0;
    recoLowPtElectronPt_.clear();
    recoLowPtElectronEta_.clear();
    recoLowPtElectronPhi_.clear();
    recoLowPtElectronVxy_.clear();
    recoLowPtElectronVz_.clear();
    recoLowPtElectronCharge_.clear();

    // Electron-positron vertices
    regreg_recoVtxVxy_.clear();
    regreg_recoVtxVz_.clear();
    regreg_recoVtxSigmaVxy_.clear();
    regreg_recoVtxReducedChi2_.clear();
    regreg_recoVtxDr_.clear();
    lowlow_recoVtxVxy_.clear();
    lowlow_recoVtxVz_.clear();
    lowlow_recoVtxSigmaVxy_.clear();
    lowlow_recoVtxReducedChi2_.clear();
    lowlow_recoVtxDr_.clear();
    lowreg_recoVtxVxy_.clear();
    lowreg_recoVtxVz_.clear();
    lowreg_recoVtxSigmaVxy_.clear();
    lowreg_recoVtxReducedChi2_.clear();
    lowreg_recoVtxDr_.clear();
}
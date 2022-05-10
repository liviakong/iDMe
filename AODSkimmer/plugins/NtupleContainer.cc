#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetTree(TTree *tree) { outT = tree; }

void NtupleContainer::CreateTreeBranches() {

    // Data or MC
    //outT->Branch("isData",&isData_);

    // Reco information
    outT->Branch("eventNum", &eventNum_);
    outT->Branch("lumiSec", &lumiSec_);
    outT->Branch("runNum", &runNum_);

    // Normal Electrons
    outT->Branch("nElectron",&nElectronDefault_);
    outT->Branch("Electron_pt",&recoElectronPt_);
    outT->Branch("Electron_eta",&recoElectronEta_);
    outT->Branch("Electron_phi",&recoElectronPhi_);
    outT->Branch("Electron_energy",&recoElectronE_);
    outT->Branch("Electron_px",&recoElectronPx_);
    outT->Branch("Electron_py",&recoElectronPy_);
    outT->Branch("Electron_pz",&recoElectronPz_);
    outT->Branch("Electron_vxy",&recoElectronVxy_);
    outT->Branch("Electron_vz",&recoElectronVz_);
    outT->Branch("Electron_dxy",&recoElectronDxy_);
    outT->Branch("Electron_dxyErr",&recoElectronDxyError_);
    outT->Branch("Electron_dz",&recoElectronDz_);
    outT->Branch("Electron_dzErr",&recoElectronDzError_);
    outT->Branch("Electron_trkChi2",&recoElectronTrkChi2_);
    outT->Branch("Electron_trkIso",&recoElectronTrkIso_);
    outT->Branch("Electron_numTrackerHits",&recoElectronTrkNumTrackerHits_);
    outT->Branch("Electron_numPixHits",&recoElectronTrkNumPixHits_);
    outT->Branch("Electron_numStripHits",&recoElectronTrkNumStripHits_);
    outT->Branch("Electron_charge",&recoElectronCharge_);
    outT->Branch("Electron_passConvVeto",&recoElectron_passConversionVeto_);

    // Low pT electrons
    outT->Branch("nLptElectron",&nElectronLowPt_);
    outT->Branch("LptElectron_pt",&recoLowPtElectronPt_);
    outT->Branch("LptElectron_eta",&recoLowPtElectronEta_);
    outT->Branch("LptElectron_phi",&recoLowPtElectronPhi_);
    outT->Branch("LptElectron_energy",&recoLowPtElectronE_);
    outT->Branch("LptElectron_px",&recoLowPtElectronPx_);
    outT->Branch("LptElectron_py",&recoLowPtElectronPy_);
    outT->Branch("LptElectron_pz",&recoLowPtElectronPz_);
    outT->Branch("LptElectron_vxy",&recoLowPtElectronVxy_);
    outT->Branch("LptElectron_vz",&recoLowPtElectronVz_);
    outT->Branch("LptElectron_dxy",&recoLowPtElectronDxy_);
    outT->Branch("LptElectron_dxyErr",&recoLowPtElectronDxyError_);
    outT->Branch("LptElectron_dz",&recoLowPtElectronDz_);
    outT->Branch("LptElectron_dzErr",&recoLowPtElectronDzError_);
    outT->Branch("LptElectron_trkChi2",&recoLowPtElectronTrkChi2_);
    outT->Branch("LptElectron_trkIso",&recoLowPtElectronTrkIso_);
    outT->Branch("LptElectron_numTrackerHits",&recoLowPtElectronTrkNumTrackerHits_);
    outT->Branch("LptElectron_numPixHits",&recoLowPtElectronTrkNumPixHits_);
    outT->Branch("LptElectron_numStripHits",&recoLowPtElectronTrkNumStripHits_);
    outT->Branch("LptElectron_charge",&recoLowPtElectronCharge_);
    outT->Branch("LptElectron_passConvVeto",&recoLowPtElectron_passConversionVeto_);

    // Photons
    outT->Branch("nPhoton",&nPhotons_);
    outT->Branch("Photon_et",&PhotonEt_);
    outT->Branch("Photon_eta",&PhotonEta_);
    outT->Branch("Photon_phi",&PhotonPhi_);

    // OOT Photons
    outT->Branch("nootPhoton",&nOOTPhotons_);
    outT->Branch("ootPhoton_et",&ootPhotonEt_);
    outT->Branch("ootPhoton_eta",&ootPhotonEta_);
    outT->Branch("ootPhoton_phi",&ootPhotonPhi_);

    // Photon conversions
    outT->Branch("nConversion",&nConversions_);
    outT->Branch("Conversion_pt",&conversionPt_);
    outT->Branch("Conversion_eta",&conversionEta_);
    outT->Branch("Conversion_phi",&conversionPhi_);
    outT->Branch("Conversion_energy",&conversionE_);
    outT->Branch("Conversion_px",&conversionPx_);
    outT->Branch("Conversion_py",&conversionPy_);
    outT->Branch("Conversion_pz",&conversionPz_);
    outT->Branch("Conversion_vxy",&conversionVxy_);
    outT->Branch("Conversion_vz",&conversionVz_);
    outT->Branch("Conversion_x",&conversionX_);
    outT->Branch("Conversion_y",&conversionY_);
    outT->Branch("Conversion_z",&conversionZ_);
    outT->Branch("Conversion_trk1_innerPt",&conversion_Trk1_innerPt_);
    outT->Branch("Conversion_trk1_innerEta",&conversion_Trk1_innerEta_);
    outT->Branch("Conversion_trk1_innerPhi",&conversion_Trk1_innerPhi_);
    outT->Branch("Conversion_trk1_outerPt",&conversion_Trk1_outerPt_);
    outT->Branch("Conversion_trk1_outerEta",&conversion_Trk1_outerEta_);
    outT->Branch("Conversion_trk1_outerPhi",&conversion_Trk1_outerPhi_);
    outT->Branch("Conversion_trk2_innerPt",&conversion_Trk2_innerPt_);
    outT->Branch("Conversion_trk2_innerEta",&conversion_Trk2_innerEta_);
    outT->Branch("Conversion_trk2_innerPhi",&conversion_Trk2_innerPhi_);
    outT->Branch("Conversion_trk2_outerPt",&conversion_Trk2_outerPt_);
    outT->Branch("Conversion_trk2_outerEta",&conversion_Trk2_outerEta_);
    outT->Branch("Conversion_trk2_outerPhi",&conversion_Trk2_outerPhi_);

    // MET
    outT->Branch("MET_PFMET_ET",&PFMET_ET_);
    outT->Branch("MET_PFMET_px",&PFMET_Px_);
    outT->Branch("MET_PFMET_py",&PFMET_Py_);
    outT->Branch("MET_PFMET_pt",&PFMET_Pt_);
    outT->Branch("MET_PFMET_phi",&PFMET_Phi_);

    outT->Branch("MET_CaloMET_ET",&CaloMET_ET_);
    outT->Branch("MET_CaloMET_px",&CaloMET_Px_);
    outT->Branch("MET_CaloMET_py",&CaloMET_Py_);
    outT->Branch("MET_CaloMET_pt",&CaloMET_Pt_);
    outT->Branch("MET_CaloMET_phi",&CaloMET_Phi_);
    
    outT->Branch("PuppiMET_PFMET_ET",&PuppiPFMET_ET_);
    outT->Branch("PuppiMET_PFMET_px",&PuppiPFMET_Px_);
    outT->Branch("PuppiMET_PFMET_py",&PuppiPFMET_Py_);
    outT->Branch("PuppiMET_PFMET_pt",&PuppiPFMET_Pt_);
    outT->Branch("PuppiMET_PFMET_phi",&PuppiPFMET_Phi_);
    
    outT->Branch("PuppiMET_CaloMET_ET",&PuppiCaloMET_ET_);
    outT->Branch("PuppiMET_CaloMET_px",&PuppiCaloMET_Px_);
    outT->Branch("PuppiMET_CaloMET_py",&PuppiCaloMET_Py_);
    outT->Branch("PuppiMET_CaloMET_pt",&PuppiCaloMET_Pt_);
    outT->Branch("PuppiMET_CaloMET_phi",&PuppiCaloMET_Phi_);


    // Electron-positron vertex branches
    outT->Branch("nRRvtx",&nEleVertex_RR_);
    outT->Branch("RRvtx_idx1", &RRvtx_idx1_);
    outT->Branch("RRvtx_idx2", &RRvtx_idx2_);
    outT->Branch("RRvtx_vxy", &RRvtx_recoVtxVxy_);
    outT->Branch("RRvtx_sigmavxy", &RRvtx_recoVtxSigmaVxy_);
    outT->Branch("RRvtx_vx",&RRvtx_recoVtxVx_);
    outT->Branch("RRvtx_vy",&RRvtx_recoVtxVy_);
    outT->Branch("RRvtx_vz",&RRvtx_recoVtxVz_);
    outT->Branch("RRvtx_reduced_chi2", &RRvtx_recoVtxReducedChi2_);
    outT->Branch("RRvtx_dR",  &RRvtx_recoVtxDr_);
    outT->Branch("RRvtx_sign",&RRvtx_recoVtxSign_);
    outT->Branch("RRvtx_pt",&RRvtx_ll_pt_);
    outT->Branch("RRvtx_eta",&RRvtx_ll_eta_);
    outT->Branch("RRvtx_phi",&RRvtx_ll_phi_);
    outT->Branch("RRvtx_energy",&RRvtx_ll_e_);
    outT->Branch("RRvtx_m",&RRvtx_ll_m_);
    outT->Branch("RRvtx_px",&RRvtx_ll_px_);
    outT->Branch("RRvtx_py",&RRvtx_ll_py_);
    outT->Branch("RRvtx_pz",&RRvtx_ll_pz_);

    outT->Branch("nLLvtx",&nEleVertex_LL_);
    outT->Branch("LLvtx_idx1", &LLvtx_idx1_);
    outT->Branch("LLvtx_idx2", &LLvtx_idx2_);
    outT->Branch("LLvtx_vxy", &LLvtx_recoVtxVxy_);
    outT->Branch("LLvtx_sigmavxy", &LLvtx_recoVtxSigmaVxy_);
    outT->Branch("LLvtx_vx",&LLvtx_recoVtxVx_);
    outT->Branch("LLvtx_vy",&LLvtx_recoVtxVy_);
    outT->Branch("LLvtx_vz",&LLvtx_recoVtxVz_);
    outT->Branch("LLvtx_reduced_chi2", &LLvtx_recoVtxReducedChi2_);
    outT->Branch("LLvtx_dR",  &LLvtx_recoVtxDr_);
    outT->Branch("LLvtx_sign",&LLvtx_recoVtxSign_);
    outT->Branch("LLvtx_pt",&LLvtx_ll_pt_);
    outT->Branch("LLvtx_eta",&LLvtx_ll_eta_);
    outT->Branch("LLvtx_phi",&LLvtx_ll_phi_);
    outT->Branch("LLvtx_energy",&LLvtx_ll_e_);
    outT->Branch("LLvtx_m",&LLvtx_ll_m_);
    outT->Branch("LLvtx_px",&LLvtx_ll_px_);
    outT->Branch("LLvtx_py",&LLvtx_ll_py_);
    outT->Branch("LLvtx_pz",&LLvtx_ll_pz_);

    outT->Branch("nLRvtx",&nEleVertex_LR_);
    outT->Branch("LRvtx_idx1", &LRvtx_idx1_);
    outT->Branch("LRvtx_idx2", &LRvtx_idx2_);
    outT->Branch("LRvtx_vxy", &LRvtx_recoVtxVxy_);
    outT->Branch("LRvtx_sigmavxy", &LRvtx_recoVtxSigmaVxy_);
    outT->Branch("LRvtx_vx",&LRvtx_recoVtxVx_);
    outT->Branch("LRvtx_vy",&LRvtx_recoVtxVy_);
    outT->Branch("LRvtx_vz",&LRvtx_recoVtxVz_);
    outT->Branch("LRvtx_reduced_chi2", &LRvtx_recoVtxReducedChi2_);
    outT->Branch("LRvtx_dR",  &LRvtx_recoVtxDr_);
    outT->Branch("LRvtx_sign",&LRvtx_recoVtxSign_);
    outT->Branch("LRvtx_pt",&LRvtx_ll_pt_);
    outT->Branch("LRvtx_eta",&LRvtx_ll_eta_);
    outT->Branch("LRvtx_phi",&LRvtx_ll_phi_);
    outT->Branch("LRvtx_energy",&LRvtx_ll_e_);
    outT->Branch("LRvtx_m",&LRvtx_ll_m_);
    outT->Branch("LRvtx_px",&LRvtx_ll_px_);
    outT->Branch("LRvtx_py",&LRvtx_ll_py_);
    outT->Branch("LRvtx_pz",&LRvtx_ll_pz_);

    outT->Branch("nRCvtx",&nEleVertex_RC_);
    outT->Branch("RCvtx_idx1", &RCvtx_idx1_);
    outT->Branch("RCvtx_idx2", &RCvtx_idx2_);
    outT->Branch("RCvtx_vxy", &RCvtx_recoVtxVxy_);
    outT->Branch("RCvtx_sigmavxy", &RCvtx_recoVtxSigmaVxy_);
    outT->Branch("RCvtx_vx",&RCvtx_recoVtxVx_);
    outT->Branch("RCvtx_vy",&RCvtx_recoVtxVy_);
    outT->Branch("RCvtx_vz",&RCvtx_recoVtxVz_);
    outT->Branch("RCvtx_reduced_chi2", &RCvtx_recoVtxReducedChi2_);
    outT->Branch("RCvtx_dR",  &RCvtx_recoVtxDr_);
    outT->Branch("RCvtx_sign",&RCvtx_recoVtxSign_);
    outT->Branch("RCvtx_pt",&RCvtx_ll_pt_);
    outT->Branch("RCvtx_eta",&RCvtx_ll_eta_);
    outT->Branch("RCvtx_phi",&RCvtx_ll_phi_);
    outT->Branch("RCvtx_energy",&RCvtx_ll_e_);
    outT->Branch("RCvtx_m",&RCvtx_ll_m_);
    outT->Branch("RCvtx_px",&RCvtx_ll_px_);
    outT->Branch("RCvtx_py",&RCvtx_ll_py_);
    outT->Branch("RCvtx_pz",&RCvtx_ll_pz_);

    outT->Branch("nLCvtx",&nEleVertex_LC_);
    outT->Branch("LCvtx_idx1", &LCvtx_idx1_);
    outT->Branch("LCvtx_idx2", &LCvtx_idx2_);
    outT->Branch("LCvtx_vxy", &LCvtx_recoVtxVxy_);
    outT->Branch("LCvtx_sigmavxy", &LCvtx_recoVtxSigmaVxy_);
    outT->Branch("LCvtx_vx",&LCvtx_recoVtxVx_);
    outT->Branch("LCvtx_vy",&LCvtx_recoVtxVy_);
    outT->Branch("LCvtx_vz",&LCvtx_recoVtxVz_);
    outT->Branch("LCvtx_reduced_chi2", &LCvtx_recoVtxReducedChi2_);
    outT->Branch("LCvtx_dR",  &LCvtx_recoVtxDr_);
    outT->Branch("LCvtx_sign",&LCvtx_recoVtxSign_);
    outT->Branch("LCvtx_pt",&LCvtx_ll_pt_);
    outT->Branch("LCvtx_eta",&LCvtx_ll_eta_);
    outT->Branch("LCvtx_phi",&LCvtx_ll_phi_);
    outT->Branch("LCvtx_energy",&LCvtx_ll_e_);
    outT->Branch("LCvtx_m",&LCvtx_ll_m_);
    outT->Branch("LCvtx_px",&LCvtx_ll_px_);
    outT->Branch("LCvtx_py",&LCvtx_ll_py_);
    outT->Branch("LCvtx_pz",&LCvtx_ll_pz_);

    // Electron candidates from isoTracks / ECAL clusters
    outT->Branch("nEleCand",&nEleCand_);
    outT->Branch("EleCand_pt",&EleCand_pt_);
    outT->Branch("EleCand_et",&EleCand_et_);
    outT->Branch("EleCand_eta",&EleCand_eta_);
    outT->Branch("EleCand_phi",&EleCand_phi_);
    outT->Branch("EleCand_dxy",&EleCand_dxy_);
    outT->Branch("EleCand_dxyErr",&EleCand_dxyError_);
    outT->Branch("EleCand_dxy_PV",&EleCand_dxy_PV_);
    outT->Branch("EleCand_dxyErr_PV",&EleCand_dxyError_PV_);
    outT->Branch("EleCand_relPFiso",&EleCand_relPFiso_);
    outT->Branch("EleCand_relTrkiso",&EleCand_relTrkiso_);
    outT->Branch("EleCand_ptDiff",&EleCand_ptDiff_);
    outT->Branch("EleCand_trkIso",&EleCand_trkIso_);
    outT->Branch("EleCand_trkChi2",&EleCand_trkChi2_);
    outT->Branch("EleCand_numTrackerHits",&EleCand_numTrackerHits_);
    outT->Branch("EleCand_numPixHits",&EleCand_numPixHits_);
    outT->Branch("EleCand_numStripHits",&EleCand_numStripHits_);

    // Displaced dileptons
    outT->Branch("ndispEE",&ndispEE_);
    outT->Branch("dispEE_maxIxy",&dispEE_maxIxy_);
    outT->Branch("dispEE_Lxy",&dispEE_Lxy_);
    outT->Branch("dispEE_Ixy",&dispEE_Ixy_);
    outT->Branch("dispEE_trackDxy",&dispEE_trackDxy_);
    outT->Branch("dispEE_trackIxy",&dispEE_trackIxy_);
    outT->Branch("dispEE_vx",&dispEE_vx_);
    outT->Branch("dispEE_vy",&dispEE_vy_);
    outT->Branch("dispEE_mass",&dispEE_mass_);
    outT->Branch("dispEE_normalizedChi2",&dispEE_normalizedChi2_);
    outT->Branch("dispEE_leadingPt",&dispEE_leadingPt_);
    outT->Branch("dispEE_subleadingPt",&dispEE_subleadingPt_);
    outT->Branch("dispEE_leadingEt",&dispEE_leadingEt_);
    outT->Branch("dispEE_subleadingEt",&dispEE_subleadingEt_);
    outT->Branch("dispEE_cosAlpha",&dispEE_cosAlpha_);
    outT->Branch("dispEE_dPhi",&dispEE_dPhi_);
    outT->Branch("dispEE_relisoA",&dispEE_relisoA_);
    outT->Branch("dispEE_relisoB",&dispEE_relisoB_);
    outT->Branch("dispEE_fromPVA",&dispEE_fromPVA_);
    outT->Branch("dispEE_fromPVB",&dispEE_fromPVB_);
    outT->Branch("dispEE_PVAssociation",&dispEE_PVAssociation_);

    // Displaced dilepton candidates
    outT->Branch("nEECand",&nEECand_);
    outT->Branch("EECand_Lxy_PV",&EECand_Lxy_PV_);
    outT->Branch("EECand_Ixy_PV",&EECand_Ixy_PV_);
    outT->Branch("EECand_Lxy_0",&EECand_Lxy_0_);
    outT->Branch("EECand_Ixy_0",&EECand_Ixy_0_);
    outT->Branch("EECand_Lxy_BS",&EECand_Lxy_BS_);
    outT->Branch("EECand_Ixy_BS",&EECand_Ixy_BS_);
    outT->Branch("EECand_trackDxy",&EECand_trackDxy_);
    outT->Branch("EECand_trackIxy",&EECand_trackIxy_);
    outT->Branch("EECand_trackDxy_PV",&EECand_trackDxy_PV_);
    outT->Branch("EECand_trackIxy_PV",&EECand_trackIxy_PV_);
    outT->Branch("EECand_trackDxy_0",&EECand_trackDxy_0_);
    outT->Branch("EECand_trackIxy_0",&EECand_trackIxy_0_);
    outT->Branch("EECand_trackDxy_BS",&EECand_trackDxy_BS_);
    outT->Branch("EECand_trackIxy_BS",&EECand_trackIxy_BS_);
    outT->Branch("EECand_vx",&EECand_vx_);
    outT->Branch("EECand_vy",&EECand_vy_);
    outT->Branch("EECand_mass",&EECand_mass_);
    outT->Branch("EECand_normalizedChi2",&EECand_normalizedChi2_);
    outT->Branch("EECand_leadingPt",&EECand_leadingPt_);
    outT->Branch("EECand_subleadingPt",&EECand_subleadingPt_);
    outT->Branch("EECand_leadingEt",&EECand_leadingEt_);
    outT->Branch("EECand_subleadingEt",&EECand_subleadingEt_);
    outT->Branch("EECand_cosAlpha",&EECand_cosAlpha_);
    outT->Branch("EECand_dR",&EECand_dR_);
    outT->Branch("EECand_dPhi",&EECand_dPhi_);
    outT->Branch("EECand_lldPhi",&EECand_lldPhi_);
    outT->Branch("EECand_relisoA",&EECand_relisoA_);
    outT->Branch("EECand_relisoB",&EECand_relisoB_);

    // Gen information
    if (!isData_) {
        outT->Branch("genWgt", &genwgt_);
        outT->Branch("nGenPart",&nGen_);
        outT->Branch("GenPart_ID", &genID_);
        outT->Branch("GenPart_charge", &genCharge_);
        outT->Branch("GenPart_pt", &genPt_);
        outT->Branch("GenPart_eta", &genEta_);
        outT->Branch("GenPart_phi", &genPhi_);
        outT->Branch("GenPart_energy", &genEn_);
        outT->Branch("GenPart_px",&genPx_);
        outT->Branch("GenPart_py",&genPy_);
        outT->Branch("GenPart_pz", &genPz_);
        outT->Branch("GenPart_vxy", &genVxy_);
        outT->Branch("GenPart_vx",&genVx_);
        outT->Branch("GenPart_vy",&genVy_);
        outT->Branch("GenPart_vz",&genVz_);
        outT->Branch("GenPart_mass", &genMass_);
        outT->Branch("nGenJet",&nGenJet_);
        outT->Branch("GenJet_pt", &genJetPt_);
        outT->Branch("GenJet_eta", &genJetEta_);
        outT->Branch("GenJet_phi", &genJetPhi_);
        outT->Branch("GenMET_pt", &genLeadMETPt_);
        outT->Branch("GenMET_phi", &genLeadMETPhi_);
        outT->Branch("GenMET_px",&genLeadMETPx_);
        outT->Branch("GenMET_py",&genLeadMETPy_);
        outT->Branch("GenMET_ET",&genLeadMETET_);
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
    genEn_.clear();
    genPx_.clear();
    genPy_.clear();
    genPz_.clear();
    genVxy_.clear();
    genVx_.clear();
    genVy_.clear();
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
    genLeadMETET_ = -9999;
    genLeadMETPx_ = -9999;
    genLeadMETPy_ = -9999;

    // Electrons
    nElectronDefault_ = 0;
    recoElectronPt_.clear();
    recoElectronEta_.clear();
    recoElectronPhi_.clear();
    recoElectronE_.clear();
    recoElectronPx_.clear();
    recoElectronPy_.clear();
    recoElectronPz_.clear();
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
    recoElectron_passConversionVeto_.clear();

    // Low pT electrons
    nElectronLowPt_ = 0;
    recoLowPtElectronPt_.clear();
    recoLowPtElectronEta_.clear();
    recoLowPtElectronPhi_.clear();
    recoLowPtElectronE_.clear();
    recoLowPtElectronPx_.clear();
    recoLowPtElectronPy_.clear();
    recoLowPtElectronPz_.clear();
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
    recoLowPtElectron_passConversionVeto_.clear();

    // Photons
    nPhotons_ = 0;
    PhotonEt_.clear();
    PhotonEta_.clear();
    PhotonPhi_.clear();

    // OOT Photons
    nOOTPhotons_ = 0;
    ootPhotonEt_.clear();
    ootPhotonEta_.clear();
    ootPhotonPhi_.clear();

    // Photon conversions
    nConversions_ = 0;
    conversionPt_.clear();
    conversionEta_.clear();
    conversionPhi_.clear();
    conversionE_.clear();
    conversionPx_.clear();
    conversionPy_.clear();
    conversionPz_.clear();
    conversionVxy_.clear();
    conversionVz_.clear();
    conversionX_.clear();
    conversionY_.clear();
    conversionZ_.clear();
    conversion_Trk1_innerPt_.clear();
    conversion_Trk1_innerEta_.clear();
    conversion_Trk1_innerPhi_.clear();
    conversion_Trk1_outerPt_.clear();
    conversion_Trk1_outerEta_.clear();
    conversion_Trk1_outerPhi_.clear();
    conversion_Trk2_innerPt_.clear();
    conversion_Trk2_innerEta_.clear();
    conversion_Trk2_innerPhi_.clear();
    conversion_Trk2_outerPt_.clear();
    conversion_Trk2_outerEta_.clear();
    conversion_Trk2_outerPhi_.clear();

    // MET 
    PFMET_ET_ = -999;
    PFMET_Px_ = -999;
    PFMET_Py_ = -999;
    PFMET_Pt_ = -999;
    PFMET_Phi_ = -999;
    
    CaloMET_ET_ = -999;
    CaloMET_Px_ = -999;
    CaloMET_Py_ = -999;
    CaloMET_Pt_ = -999;
    CaloMET_Phi_ = -999;

    PuppiPFMET_ET_ = -999;
    PuppiPFMET_Px_ = -999;
    PuppiPFMET_Py_ = -999;
    PuppiPFMET_Pt_ = -999;
    PuppiPFMET_Phi_ = -999;
    
    PuppiCaloMET_ET_ = -999;
    PuppiCaloMET_Px_ = -999;
    PuppiCaloMET_Py_ = -999;
    PuppiCaloMET_Pt_ = -999;
    PuppiCaloMET_Phi_ = -999;

    // Electron-positron vertices
    nEleVertex_RR_ = 0;
    RRvtx_idx1_.clear();
    RRvtx_idx2_.clear();
    RRvtx_recoVtxVxy_.clear();
    RRvtx_recoVtxVz_.clear();
    RRvtx_recoVtxSigmaVxy_.clear();
    RRvtx_recoVtxVx_.clear();
    RRvtx_recoVtxVy_.clear();
    RRvtx_recoVtxVz_.clear();
    RRvtx_recoVtxReducedChi2_.clear();
    RRvtx_recoVtxDr_.clear();
    RRvtx_recoVtxSign_.clear();
    RRvtx_ll_pt_.clear();
    RRvtx_ll_eta_.clear();
    RRvtx_ll_phi_.clear();
    RRvtx_ll_e_.clear();
    RRvtx_ll_m_.clear();
    RRvtx_ll_px_.clear();
    RRvtx_ll_py_.clear();
    RRvtx_ll_pz_.clear();

    nEleVertex_LL_ = 0;
    LLvtx_idx1_.clear();
    LLvtx_idx2_.clear();
    LLvtx_recoVtxVxy_.clear();
    LLvtx_recoVtxVz_.clear();
    LLvtx_recoVtxSigmaVxy_.clear();
    LLvtx_recoVtxVx_.clear();
    LLvtx_recoVtxVy_.clear();
    LLvtx_recoVtxVz_.clear();
    LLvtx_recoVtxReducedChi2_.clear();
    LLvtx_recoVtxDr_.clear();
    LLvtx_recoVtxSign_.clear();
    LLvtx_ll_pt_.clear();
    LLvtx_ll_eta_.clear();
    LLvtx_ll_phi_.clear();
    LLvtx_ll_e_.clear();
    LLvtx_ll_m_.clear();
    LLvtx_ll_px_.clear();
    LLvtx_ll_py_.clear();
    LLvtx_ll_pz_.clear();

    nEleVertex_LR_ = 0;
    LRvtx_idx1_.clear();
    LRvtx_idx2_.clear();
    LRvtx_recoVtxVxy_.clear();
    LRvtx_recoVtxVz_.clear();
    LRvtx_recoVtxSigmaVxy_.clear();
    LRvtx_recoVtxVx_.clear();
    LRvtx_recoVtxVy_.clear();
    LRvtx_recoVtxVz_.clear();
    LRvtx_recoVtxReducedChi2_.clear();
    LRvtx_recoVtxDr_.clear();
    LRvtx_recoVtxSign_.clear();
    LRvtx_ll_pt_.clear();
    LRvtx_ll_eta_.clear();
    LRvtx_ll_phi_.clear();
    LRvtx_ll_e_.clear();
    LRvtx_ll_m_.clear();
    LRvtx_ll_px_.clear();
    LRvtx_ll_py_.clear();
    LRvtx_ll_pz_.clear();

    nEleVertex_RC_ = 0;
    RCvtx_idx1_.clear();
    RCvtx_idx2_.clear();
    RCvtx_recoVtxVxy_.clear();
    RCvtx_recoVtxVz_.clear();
    RCvtx_recoVtxSigmaVxy_.clear();
    RCvtx_recoVtxVx_.clear();
    RCvtx_recoVtxVy_.clear();
    RCvtx_recoVtxVz_.clear();
    RCvtx_recoVtxReducedChi2_.clear();
    RCvtx_recoVtxDr_.clear();
    RCvtx_recoVtxSign_.clear();
    RCvtx_ll_pt_.clear();
    RCvtx_ll_eta_.clear();
    RCvtx_ll_phi_.clear();
    RCvtx_ll_e_.clear();
    RCvtx_ll_m_.clear();
    RCvtx_ll_px_.clear();
    RCvtx_ll_py_.clear();
    RCvtx_ll_pz_.clear();

    nEleVertex_LC_ = 0;
    LCvtx_idx1_.clear();
    LCvtx_idx2_.clear();
    LCvtx_recoVtxVxy_.clear();
    LCvtx_recoVtxVz_.clear();
    LCvtx_recoVtxSigmaVxy_.clear();
    LCvtx_recoVtxVx_.clear();
    LCvtx_recoVtxVy_.clear();
    LCvtx_recoVtxVz_.clear();
    LCvtx_recoVtxReducedChi2_.clear();
    LCvtx_recoVtxDr_.clear();
    LCvtx_recoVtxSign_.clear();
    LCvtx_ll_pt_.clear();
    LCvtx_ll_eta_.clear();
    LCvtx_ll_phi_.clear();
    LCvtx_ll_e_.clear();
    LCvtx_ll_m_.clear();
    LCvtx_ll_px_.clear();
    LCvtx_ll_py_.clear();
    LCvtx_ll_pz_.clear();

    // Electron candidates from isotracks / ECAL clusters in displaced dielectron algo
    nEleCand_ = 0;
    EleCand_pt_.clear();
    EleCand_et_.clear();
    EleCand_eta_.clear();
    EleCand_phi_.clear();
    EleCand_dxy_.clear();
    EleCand_dxyError_.clear();
    EleCand_dxy_PV_.clear();
    EleCand_dxyError_PV_.clear();
    EleCand_relPFiso_.clear();
    EleCand_relTrkiso_.clear();
    EleCand_ptDiff_.clear();
    EleCand_trkIso_.clear();
    EleCand_trkChi2_.clear();
    EleCand_numTrackerHits_.clear();
    EleCand_numPixHits_.clear();
    EleCand_numStripHits_.clear();

    // Displaced dileptons fromd dedicated algorithm
    ndispEE_ = 0;
    dispEE_maxIxy_.clear();
    dispEE_Lxy_.clear();
    dispEE_Ixy_.clear();
    dispEE_trackDxy_.clear();
    dispEE_trackIxy_.clear();
    dispEE_vx_.clear();
    dispEE_vy_.clear();
    dispEE_mass_.clear();
    dispEE_normalizedChi2_.clear();
    dispEE_leadingPt_.clear();
    dispEE_subleadingPt_.clear();
    dispEE_leadingEt_.clear();
    dispEE_subleadingEt_.clear();
    dispEE_cosAlpha_.clear();
    dispEE_dPhi_.clear();
    dispEE_relisoA_.clear();
    dispEE_relisoB_.clear();
    dispEE_fromPVA_.clear();
    dispEE_fromPVB_.clear();
    dispEE_PVAssociation_.clear();

    // Displaced dilepton candidates
    nEECand_ = 0;
    EECand_Lxy_PV_.clear();
    EECand_Ixy_PV_.clear();
    EECand_Lxy_0_.clear();
    EECand_Ixy_0_.clear();
    EECand_Lxy_BS_.clear();
    EECand_Ixy_BS_.clear();
    EECand_trackDxy_.clear();
    EECand_trackIxy_.clear();
    EECand_trackDxy_PV_.clear();
    EECand_trackIxy_PV_.clear();
    EECand_trackDxy_0_.clear();
    EECand_trackIxy_0_.clear();
    EECand_trackDxy_BS_.clear();
    EECand_trackIxy_BS_.clear();
    EECand_vx_.clear();
    EECand_vy_.clear();
    EECand_mass_.clear();
    EECand_normalizedChi2_.clear();
    EECand_leadingPt_.clear();
    EECand_subleadingPt_.clear();
    EECand_leadingEt_.clear();
    EECand_subleadingEt_.clear();
    EECand_cosAlpha_.clear();
    EECand_dR_.clear();
    EECand_dPhi_.clear();
    EECand_lldPhi_.clear();
    EECand_relisoA_.clear();
    EECand_relisoB_.clear();
}
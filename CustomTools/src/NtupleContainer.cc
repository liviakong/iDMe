#include "iDMe/CustomTools/interface/NtupleContainer.hh"

NtupleContainer::NtupleContainer() {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetTree(TTree *tree) { outT = tree; }

void NtupleContainer::CreateTreeBranches() {

    // Data or MC
    //outT->Branch("isData",&isData_);

    // Reco information
    outT->Branch("trigFired",&fired_);
    outT->Branch("trigFired16",&fired16_);
    outT->Branch("trigFired17",&fired17_);
    outT->Branch("trigFired18",&fired18_);
    outT->Branch("eventNum", &eventNum_);
    outT->Branch("lumiSec", &lumiSec_);
    outT->Branch("runNum", &runNum_);

    // MET Filters
    outT->Branch("METFiltersFailBits",&METFiltersFailBits_);

    // Normal Electrons
    outT->Branch("nElectron",&nElectronDefault_);
    outT->Branch("Electron_matched",&recoElectronMatched_);
    outT->Branch("Electron_matchType",&recoElectronMatchType_);
    outT->Branch("Electron_matchdR",&recoElectronMatchdR_);
    outT->Branch("Electron_pt",&recoElectronPt_);
    outT->Branch("Electron_eta",&recoElectronEta_);
    outT->Branch("Electron_etaErr",&recoElectronEtaError_);
    outT->Branch("Electron_phi",&recoElectronPhi_);
    outT->Branch("Electron_phiErr",&recoElectronPhiError_);
    outT->Branch("Electron_ID",&recoElectronID_);
    outT->Branch("Electron_angRes",&recoElectronAngularRes_);
    outT->Branch("Electron_e",&recoElectronE_);
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
    outT->Branch("Electron_trkRelIso",&recoElectronTrkRelIso_);
    outT->Branch("Electron_trkProb",&recoElectronTrkProb_);
    outT->Branch("Electron_numTrackerHits",&recoElectronTrkNumTrackerHits_);
    outT->Branch("Electron_numPixHits",&recoElectronTrkNumPixHits_);
    outT->Branch("Electron_numStripHits",&recoElectronTrkNumStripHits_);
    outT->Branch("Electron_charge",&recoElectronCharge_);

    // Low pT electrons
    outT->Branch("nLptElectron",&nElectronLowPt_);
    outT->Branch("LptElectron_matched",&recoLowPtElectronMatched_);
    outT->Branch("LptElectron_matchType",&recoLowPtElectronMatchType_);
    outT->Branch("LptElectron_matchdR",&recoLowPtElectronMatchdR_);
    outT->Branch("LptElectron_pt",&recoLowPtElectronPt_);
    outT->Branch("LptElectron_eta",&recoLowPtElectronEta_);
    outT->Branch("LptElectron_etaErr",&recoLowPtElectronEtaError_);
    outT->Branch("LptElectron_phi",&recoLowPtElectronPhi_);
    outT->Branch("LptElectron_phiErr",&recoLowPtElectronPhiError_);
    outT->Branch("LptElectron_ID",&recoLowPtElectronID_);
    outT->Branch("LptElectron_angRes",&recoLowPtElectronAngularRes_);
    outT->Branch("LptElectron_e",&recoLowPtElectronE_);
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
    outT->Branch("LptElectron_trkRelIso",&recoLowPtElectronTrkRelIso_);
    outT->Branch("LptElectron_trkProb",&recoLowPtElectronTrkProb_);
    outT->Branch("LptElectron_numTrackerHits",&recoLowPtElectronTrkNumTrackerHits_);
    outT->Branch("LptElectron_numPixHits",&recoLowPtElectronTrkNumPixHits_);
    outT->Branch("LptElectron_numStripHits",&recoLowPtElectronTrkNumStripHits_);
    outT->Branch("LptElectron_charge",&recoLowPtElectronCharge_);
    outT->Branch("LptElectron_minDRtoReg",&recoLowPtElectronMinDrToReg_);

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

    // Jets
    outT->Branch("nPFJetAll",&PFNJet_);
    outT->Branch("nPFJetPassID",&PFNPassIDJet_);
    outT->Branch("nPFJet",&PFNHighPtJet_);
    outT->Branch("PFJet_pt",&PFJetPt_);
    outT->Branch("PFJet_eta",&PFJetEta_);
    outT->Branch("PFJet_phi",&PFJetPhi_);
    outT->Branch("PFJet_corrPt",&PFJetCorrectedPt_);
    outT->Branch("PFJet_corrEta",&PFJetCorrectedEta_);
    outT->Branch("PFJet_corrPhi",&PFJetCorrectedPhi_);
    outT->Branch("PFJet_bTag",&PFJetCorrectedBTag_);
    outT->Branch("PFJet_CHEF",&PFJetCorrectedCHEF_);
    outT->Branch("PFJet_NHEF",&PFJetCorrectedNHEF_);
    outT->Branch("PFJet_CEEF",&PFJetCorrectedCEEF_);
    outT->Branch("PFJet_NEEF",&PFJetCorrectedNEEF_);
    outT->Branch("PFJet_corrNumDaughters",&PFJetCorrectedNumDaughters_);
    outT->Branch("PFJet_corrCHM",&PFJetCorrectedChargedMultiplicity_);
    outT->Branch("PFJet_corrJESUp_pt",&PFJetCorrectedJESUpPt_);
    outT->Branch("PFJet_corrJESUp_eta",&PFJetCorrectedJESUpEta_);
    outT->Branch("PFJet_corrJESUp_phi",&PFJetCorrectedJESUpPhi_);
    outT->Branch("PFJet_corrJESDown_pt",&PFJetCorrectedJESDownPt_);
    outT->Branch("PFJet_corrJESDown_eta",&PFJetCorrectedJESDownEta_);
    outT->Branch("PFJet_corrJESDown_phi",&PFJetCorrectedJESDownPhi_);
    outT->Branch("PFJet_corrJERUp_pt",&PFJetCorrectedJERUpPt_);
    outT->Branch("PFJet_corrJERUp_eta",&PFJetCorrectedJERUpEta_);
    outT->Branch("PFJet_corrJERUp_phi",&PFJetCorrectedJERUpPhi_);
    outT->Branch("PFJet_corrJERDown_pt",&PFJetCorrectedJERDownPt_);
    outT->Branch("PFJet_corrJERDown_eta",&PFJetCorrectedJERDownEta_);
    outT->Branch("PFJet_corrJERDown_phi",&PFJetCorrectedJERDownPhi_);
    outT->Branch("HEM_flag",&PFHEMFlag_);

    // MET
    outT->Branch("PFMET_ET",&PFMET_ET_);
    outT->Branch("PFMET_px",&PFMET_Px_);
    outT->Branch("PFMET_py",&PFMET_Py_);
    outT->Branch("PFMET_pt",&PFMET_Pt_);
    outT->Branch("PFMET_phi",&PFMET_Phi_);
    outT->Branch("PFMET_smearingPt",&PFMETSmearingOnlyPt_);
    outT->Branch("PFMET_smearingPhi",&PFMETSmearingOnlyPhi_);
    outT->Branch("PFMET_correctedPt",&PFMETCorrectedPt_);
    outT->Branch("PFMET_correctedPhi",&PFMETCorrectedPhi_);
    outT->Branch("PFMET_deltaPx",&PFMETEEDeltaPx_);
    outT->Branch("PFMET_deltaPy",&PFMETEEDeltaPy_);
    outT->Branch("PFMET_JESUpPt",&PFMETJESUpPt_);
    outT->Branch("PFMET_JESUpPhi",&PFMETJESUpPhi_);
    outT->Branch("PFMET_JESDownPt",&PFMETJESDownPt_);
    outT->Branch("PFMET_JESDownPhi",&PFMETJESDownPhi_);
    outT->Branch("PFMET_JERUpPt",&PFMETJERUpPt_);
    outT->Branch("PFMET_JERUpPhi",&PFMETJERUpPhi_);
    outT->Branch("PFMET_JERDownPt",&PFMETJERDownPt_);
    outT->Branch("PFMET_JERDownPhi",&PFMETJERDownPhi_);
    outT->Branch("PFMET_muonEtFrac",&PFMETMuonEtFraction_);

    outT->Branch("CaloMET_ET",&CaloMET_ET_);
    outT->Branch("CaloMET_px",&CaloMET_Px_);
    outT->Branch("CaloMET_py",&CaloMET_Py_);
    outT->Branch("CaloMET_pt",&CaloMET_Pt_);
    outT->Branch("CaloMET_phi",&CaloMET_Phi_);
    
    outT->Branch("PuppiMETPF_ET",&PuppiPFMET_ET_);
    outT->Branch("PuppiMETPF_px",&PuppiPFMET_Px_);
    outT->Branch("PuppiMETPF_py",&PuppiPFMET_Py_);
    outT->Branch("PuppiMETPF_pt",&PuppiPFMET_Pt_);
    outT->Branch("PuppiMETPF_phi",&PuppiPFMET_Phi_);
    
    outT->Branch("PuppiMETCalo_ET",&PuppiCaloMET_ET_);
    outT->Branch("PuppiMETCalo_px",&PuppiCaloMET_Px_);
    outT->Branch("PuppiMETCalo_py",&PuppiCaloMET_Py_);
    outT->Branch("PuppiMETCalo_pt",&PuppiCaloMET_Pt_);
    outT->Branch("PuppiMETCalo_phi",&PuppiCaloMET_Phi_);

    // Pileup density
    outT->Branch("rho",&rho_);

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
    outT->Branch("RRvtx_prob", &RRvtx_prob_);
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
    outT->Branch("LLvtx_prob", &LLvtx_prob_);
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
    outT->Branch("LRvtx_prob", &LRvtx_prob_);
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

    // Gen information
    if (!isData_) {
        outT->Branch("genWgt", &genwgt_);
        outT->Branch("genPU_obs", &genpuobs_);
        outT->Branch("genPU_true", &genputrue_);
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
    if (!isData_ && isSignal_) {
        outT->Branch("nGenPart",&nGen_);
        outT->Branch("GenPart_ID", &genID_);
        outT->Branch("GenPart_charge", &genCharge_);
        outT->Branch("GenPart_pt", &genPt_);
        outT->Branch("GenPart_eta", &genEta_);
        outT->Branch("GenPart_phi", &genPhi_);
        outT->Branch("GenPart_e", &genEn_);
        outT->Branch("GenPart_px",&genPx_);
        outT->Branch("GenPart_py",&genPy_);
        outT->Branch("GenPart_pz", &genPz_);
        outT->Branch("GenPart_vxy", &genVxy_);
        outT->Branch("GenPart_vx",&genVx_);
        outT->Branch("GenPart_vy",&genVy_);
        outT->Branch("GenPart_vz",&genVz_);
        outT->Branch("GenPart_mass", &genMass_);

        outT->Branch("GenEle_charge",&genEleCharge_);
        outT->Branch("GenEle_motherID",&genEleMotherID_);
        outT->Branch("GenEle_pt",&genElePt_);
        outT->Branch("GenEle_eta",&genEleEta_);
        outT->Branch("GenEle_phi",&genElePhi_);
        outT->Branch("GenEle_energy",&genEleEn_);
        outT->Branch("GenEle_px",&genElePx_);
        outT->Branch("GenEle_py",&genElePy_);
        outT->Branch("GenEle_pz",&genElePz_);
        outT->Branch("GenEle_vxy",&genEleVxy_);
        outT->Branch("GenEle_vx",&genEleVz_);
        outT->Branch("GenEle_vy",&genEleVx_);
        outT->Branch("GenEle_vz",&genEleVy_);
        outT->Branch("GenEleClosest_typ",&genEleClosestType_); // need to misspell type so it doesn't cause python issues
        outT->Branch("GenEleClosest_ind",&genEleClosestInd_);
        outT->Branch("GenEleClosest_dr",&genEleClosestDr_);
        outT->Branch("GenEleClosestReg_ind",&genEleClosestInd_reg_);
        outT->Branch("GenEleClosestReg_dr",&genEleClosestDr_reg_);
        outT->Branch("GenEleClosestLpt_ind",&genEleClosestInd_lpt_);
        outT->Branch("GenEleClosestLpt_dr",&genEleClosestDr_lpt_);

        outT->Branch("GenPos_charge",&genPosCharge_);
        outT->Branch("GenPos_motherID",&genPosMotherID_);
        outT->Branch("GenPos_pt",&genPosPt_);
        outT->Branch("GenPos_eta",&genPosEta_);
        outT->Branch("GenPos_phi",&genPosPhi_);
        outT->Branch("GenPos_energy",&genPosEn_);
        outT->Branch("GenPos_px",&genPosPx_);
        outT->Branch("GenPos_py",&genPosPy_);
        outT->Branch("GenPos_pz",&genPosPz_);
        outT->Branch("GenPos_vxy",&genPosVxy_);
        outT->Branch("GenPos_vx",&genPosVz_);
        outT->Branch("GenPos_vy",&genPosVx_);
        outT->Branch("GenPos_vz",&genPosVy_);
        outT->Branch("GenPosClosest_typ",&genPosClosestType_); // need to misspell type so it doesn't cause python issues
        outT->Branch("GenPosClosest_ind",&genPosClosestInd_);
        outT->Branch("GenPosClosest_dr",&genPosClosestDr_);
        outT->Branch("GenPosClosestReg_ind",&genPosClosestInd_reg_);
        outT->Branch("GenPosClosestReg_dr",&genPosClosestDr_reg_);
        outT->Branch("GenPosClosestLpt_ind",&genPosClosestInd_lpt_);
        outT->Branch("GenPosClosestLpt_dr",&genPosClosestDr_lpt_);

        // Gen Electron + Positron info
        outT->Branch("genEE_pt",&genEEPt_);
        outT->Branch("genEE_eta",&genEEEta_);
        outT->Branch("genEE_phi",&genEEPhi_);
        outT->Branch("genEE_energy",&genEEEn_);
        outT->Branch("genEE_mass",&genEEMass_);
        outT->Branch("genEE_dr",&genEEdR_);
    }

}

void NtupleContainer::ClearTreeBranches() {
    // Reset trigger
    fired_ = 0;
    fired16_ = 0;
    fired17_ = 0;
    fired18_ = 0;

    // MET Filters
    METFiltersFailBits_ = 0;

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

    // Gen Electron & Positron from iDM signal
    genEleCharge_ = 0;
    genEleMotherID_ = 0;
    genElePt_ = -999;
    genEleEta_ = -999;
    genElePhi_ = -999;
    genEleEn_ = -999;
    genElePx_ = -999;
    genElePy_ = -999;
    genElePz_ = -999;
    genEleVxy_ = -999;
    genEleVz_ = -999;
    genEleVx_ = -999;
    genEleVy_ = -999;
    genEleClosestType_ = -1;
    genEleClosestInd_ = -1;
    genEleClosestDr_ = 999.0;
    genEleClosestInd_reg_ = -1;
    genEleClosestDr_reg_ = 999.0;
    genEleClosestInd_lpt_ = -1;
    genEleClosestDr_lpt_ = 999.0;

    genPosCharge_ = 0;
    genPosMotherID_ = 0;
    genPosPt_ = -999;
    genPosEta_ = -999;
    genPosPhi_ = -999;
    genPosEn_ = -999;
    genPosPx_ = -999;
    genPosPy_ = -999;
    genPosPz_ = -999;
    genPosVxy_ = -999;
    genPosVz_ = -999;
    genPosVx_ = -999;
    genPosVy_ = -999;
    genPosClosestType_ = -1;
    genPosClosestInd_ = -1;
    genPosClosestDr_ = 999.0;
    genPosClosestInd_reg_ = -1;
    genPosClosestDr_reg_ = 999.0;
    genPosClosestInd_lpt_ = -1;
    genPosClosestDr_lpt_ = 999.0;

    // Gen Electron + Positron info
    genEEPt_ = -999;
    genEEEta_ = -999;
    genEEPhi_ = -999;
    genEEEn_ = -999;
    genEEMass_ = -999;
    genEEdR_ = -999;
    
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
    recoElectronEtaError_.clear();
    recoElectronPhi_.clear();
    recoElectronPhiError_.clear();
    recoElectronID_.clear();
    recoElectronAngularRes_.clear();
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
    recoElectronTrkRelIso_.clear();
    recoElectronTrkProb_.clear();
    recoElectronTrkNumTrackerHits_.clear();
    recoElectronTrkNumPixHits_.clear();
    recoElectronTrkNumStripHits_.clear();
    recoElectronCharge_.clear();

    // Low pT electrons
    nElectronLowPt_ = 0;
    recoLowPtElectronPt_.clear();
    recoLowPtElectronEta_.clear();
    recoLowPtElectronEtaError_.clear();
    recoLowPtElectronPhi_.clear();
    recoLowPtElectronPhiError_.clear();
    recoLowPtElectronID_.clear();
    recoLowPtElectronAngularRes_.clear();
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
    recoLowPtElectronTrkRelIso_.clear();
    recoLowPtElectronTrkProb_.clear();
    recoLowPtElectronTrkNumTrackerHits_.clear();
    recoLowPtElectronTrkNumPixHits_.clear();
    recoLowPtElectronTrkNumStripHits_.clear();
    recoLowPtElectronCharge_.clear();
    recoLowPtElectronMinDrToReg_.clear();

    // Gen weight and pileup
    genwgt_ = -9999;
    genpuobs_ = -9999;
    genputrue_ = -9999;

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

    // Jets
    PFNJet_ = 0;
    PFNPassIDJet_ = 0;
    PFNHighPtJet_ = 0;
    PFJetPt_.clear();
    PFJetEta_.clear();
    PFJetPhi_.clear();
    PFJetCorrectedPt_.clear();
    PFJetCorrectedEta_.clear();
    PFJetCorrectedPhi_.clear();
    PFJetCorrectedBTag_.clear();
    PFJetCorrectedCHEF_.clear();
    PFJetCorrectedNHEF_.clear();
    PFJetCorrectedCEEF_.clear();
    PFJetCorrectedNEEF_.clear();
    PFJetCorrectedNumDaughters_.clear();
    PFJetCorrectedChargedMultiplicity_.clear();
    PFJetCorrectedJESUpPt_.clear();
    PFJetCorrectedJESUpEta_.clear();
    PFJetCorrectedJESUpPhi_.clear();
    PFJetCorrectedJESDownPt_.clear();
    PFJetCorrectedJESDownEta_.clear();
    PFJetCorrectedJESDownPhi_.clear();
    PFJetCorrectedJERUpPt_.clear();
    PFJetCorrectedJERUpEta_.clear();
    PFJetCorrectedJERUpPhi_.clear();
    PFJetCorrectedJERDownPt_.clear();
    PFJetCorrectedJERDownEta_.clear();
    PFJetCorrectedJERDownPhi_.clear();
    PFHEMFlag_ = false;

    // MET 
    PFMET_ET_ = -999;
    PFMET_Px_ = -999;
    PFMET_Py_ = -999;
    PFMET_Pt_ = -999;
    PFMET_Phi_ = -999;
    PFMETSmearingOnlyPt_ = -9999;
    PFMETSmearingOnlyPhi_ = -9999;
    PFMETCorrectedPt_ = -9999;
    PFMETCorrectedPhi_ = -9999;
    PFMETEEDeltaPx_ = 0.0;
    PFMETEEDeltaPy_ = 0.0;
    PFMETJESUpPt_ = -9999;
    PFMETJESUpPhi_ = -9999;
    PFMETJESDownPt_ = -9999;
    PFMETJESDownPhi_ = -9999;
    PFMETJERUpPt_ = -9999;
    PFMETJERUpPhi_ = -9999;
    PFMETJERDownPt_ = -9999;
    PFMETJERDownPhi_ = -9999;
    PFMETMuonEtFraction_ = -9999;
    
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

    // Pileup density
    rho_ = -9999;

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
    RRvtx_prob_.clear();
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
    LLvtx_prob_.clear();
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
    LRvtx_prob_.clear();
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
}
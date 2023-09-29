#include "iDMe/CustomTools/interface/NtupleContainerV2.hh"

NtupleContainerV2::NtupleContainerV2() {}

NtupleContainerV2::~NtupleContainerV2() {}

void NtupleContainerV2::SetTree(TTree *tree) { outT = tree; }

void NtupleContainerV2::CreateTreeBranches() {

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
    outT->Branch("Electron_pt",&recoElectronPt_);
    outT->Branch("Electron_eta",&recoElectronEta_);
    outT->Branch("Electron_etaErr",&recoElectronEtaError_);
    outT->Branch("Electron_phi",&recoElectronPhi_);
    outT->Branch("Electron_phiErr",&recoElectronPhiError_);
    outT->Branch("Electron_IDcutVeto",&recoElectronID_cutVeto_);
    outT->Branch("Electron_IDcutLoose",&recoElectronID_cutLoose_);
    outT->Branch("Electron_IDcutMed",&recoElectronID_cutMed_);
    outT->Branch("Electron_IDcutTight",&recoElectronID_cutTight_);
    outT->Branch("Electron_IDmvaIso90",&recoElectronID_mvaIso90_);
    outT->Branch("Electron_IDmvaIso80",&recoElectronID_mvaIso80_);
    outT->Branch("Electron_IDmvaIsoLoose",&recoElectronID_mvaIsoLoose_);
    outT->Branch("Electron_IDmva90",&recoElectronID_mva90_);
    outT->Branch("Electron_IDmva80",&recoElectronID_mva80_);
    outT->Branch("Electron_IDmvaLoose",&recoElectronID_mvaLoose_);
    outT->Branch("Electron_angRes",&recoElectronAngularRes_);
    outT->Branch("Electron_e",&recoElectronE_);
    outT->Branch("Electron_vxy",&recoElectronVxy_);
    outT->Branch("Electron_vz",&recoElectronVz_);
    outT->Branch("Electron_dxy",&recoElectronDxy_);
    outT->Branch("Electron_dxyErr",&recoElectronDxyError_);
    outT->Branch("Electron_dz",&recoElectronDz_);
    outT->Branch("Electron_dzErr",&recoElectronDzError_);
    outT->Branch("Electron_trkChi2",&recoElectronTrkChi2_);
    outT->Branch("Electron_trkIso",&recoElectronTrkIso_);
    outT->Branch("Electron_trkRelIso",&recoElectronTrkRelIso_);
    outT->Branch("Electron_calIso",&recoElectronCaloIso_);
    outT->Branch("Electron_calRelIso",&recoElectronCaloRelIso_);
    outT->Branch("Electron_PFIso4",&recoElectronPFIso_dR4_);
    outT->Branch("Electron_PFRelIso4",&recoElectronPFRelIso_dR4_);
    outT->Branch("Electron_PFIso3",&recoElectronPFIso_dR3_);
    outT->Branch("Electron_PFRelIso3",&recoElectronPFRelIso_dR3_);
    outT->Branch("Electron_PFIso8",&recoElectronPFIso_dR8_);
    outT->Branch("Electron_PFRelIso8",&recoElectronPFRelIso_dR8_);
    outT->Branch("Electron_PFIso",&recoElectronPFIso_);
    outT->Branch("Electron_PFRelIso",&recoElectronPFRelIso_);
    outT->Branch("Electron_chadIso",&recoElectronChadIso_);
    outT->Branch("Electron_nhadIso",&recoElectronNhadIso_);
    outT->Branch("Electron_phoIso",&recoElectronPhoIso_);
    outT->Branch("Electron_rhoEA",&recoElectronRhoEA_);
    outT->Branch("Electron_trkProb",&recoElectronTrkProb_);
    outT->Branch("Electron_numTrackerHits",&recoElectronTrkNumTrackerHits_);
    outT->Branch("Electron_numPixHits",&recoElectronTrkNumPixHits_);
    outT->Branch("Electron_numStripHits",&recoElectronTrkNumStripHits_);
    outT->Branch("Electron_charge",&recoElectronCharge_);

    // Low pT electrons
    outT->Branch("nLptElectron",&nElectronLowPt_);
    outT->Branch("LptElectron_pt",&recoLowPtElectronPt_);
    outT->Branch("LptElectron_eta",&recoLowPtElectronEta_);
    outT->Branch("LptElectron_etaErr",&recoLowPtElectronEtaError_);
    outT->Branch("LptElectron_phi",&recoLowPtElectronPhi_);
    outT->Branch("LptElectron_phiErr",&recoLowPtElectronPhiError_);
    outT->Branch("LptElectron_ID",&recoLowPtElectronID_);
    outT->Branch("LptElectron_angRes",&recoLowPtElectronAngularRes_);
    outT->Branch("LptElectron_e",&recoLowPtElectronE_);
    outT->Branch("LptElectron_vxy",&recoLowPtElectronVxy_);
    outT->Branch("LptElectron_vz",&recoLowPtElectronVz_);
    outT->Branch("LptElectron_dxy",&recoLowPtElectronDxy_);
    outT->Branch("LptElectron_dxyErr",&recoLowPtElectronDxyError_);
    outT->Branch("LptElectron_dz",&recoLowPtElectronDz_);
    outT->Branch("LptElectron_dzErr",&recoLowPtElectronDzError_);
    outT->Branch("LptElectron_trkChi2",&recoLowPtElectronTrkChi2_);
    outT->Branch("LptElectron_trkIso",&recoLowPtElectronTrkIso_);
    outT->Branch("LptElectron_trkRelIso",&recoLowPtElectronTrkRelIso_);
    outT->Branch("LptElectron_calIso",&recoLowPtElectronCaloIso_);
    outT->Branch("LptElectron_calRelIso",&recoLowPtElectronCaloRelIso_);
    outT->Branch("LptElectron_PFIso4",&recoLowPtElectronPFIso_dR4_);
    outT->Branch("LptElectron_PFRelIso4",&recoLowPtElectronPFRelIso_dR4_);
    outT->Branch("LptElectron_PFIso3",&recoLowPtElectronPFIso_dR3_);
    outT->Branch("LptElectron_PFRelIso3",&recoLowPtElectronPFRelIso_dR3_);
    outT->Branch("LptElectron_PFIso8",&recoLowPtElectronPFIso_dR8_);
    outT->Branch("LptElectron_PFIso",&recoLowPtElectronPFIso_);
    outT->Branch("LptElectron_PFRelIso",&recoLowPtElectronPFRelIso_);
    outT->Branch("LptElectron_PFRelIso8",&recoLowPtElectronPFRelIso_dR8_);
    outT->Branch("LptElectron_chadIso",&recoLowPtElectronChadIso_);
    outT->Branch("LptElectron_nhadIso",&recoLowPtElectronNhadIso_);
    outT->Branch("LptElectron_phoIso",&recoLowPtElectronPhoIso_);
    outT->Branch("LptElectron_rhoEA",&recoLowPtElectronRhoEA_);
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
    outT->Branch("nPFJetAll",&PFNJetAll_);
    outT->Branch("nPFJet",&PFNJet_);
    outT->Branch("PFJet_pt",&PFJetPt_);
    outT->Branch("PFJet_eta",&PFJetEta_);
    outT->Branch("PFJet_phi",&PFJetPhi_);
    outT->Branch("PFJet_bTag",&PFJetBTag_);
    outT->Branch("PFJet_METdPhi",&PFJetMETdPhi_);
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
    outT->Branch("PFMET_pt",&PFMET_Pt_);
    outT->Branch("PFMET_phi",&PFMET_Phi_);
    outT->Branch("PFMET_JESUpPt",&PFMETJESUpPt_);
    outT->Branch("PFMET_JESUpPhi",&PFMETJESUpPhi_);
    outT->Branch("PFMET_JESDownPt",&PFMETJESDownPt_);
    outT->Branch("PFMET_JESDownPhi",&PFMETJESDownPhi_);
    outT->Branch("PFMET_JERUpPt",&PFMETJERUpPt_);
    outT->Branch("PFMET_JERUpPhi",&PFMETJERUpPhi_);
    outT->Branch("PFMET_JERDownPt",&PFMETJERDownPt_);
    outT->Branch("PFMET_JERDownPhi",&PFMETJERDownPhi_);
    outT->Branch("PFMET_JetResUpSmearPt",&PFMETJetResUpSmearPt_);
    outT->Branch("PFMET_JetResUpSmearPhi",&PFMETJetResUpSmearPhi_);
    outT->Branch("PFMET_JetResDownSmearPt",&PFMETJetResDownSmearPt_);
    outT->Branch("PFMET_JetResDownSmearPhi",&PFMETJetResDownSmearPhi_);

    outT->Branch("CaloMET_ET",&CaloMET_ET_);
    outT->Branch("CaloMET_pt",&CaloMET_Pt_);
    outT->Branch("CaloMET_phi",&CaloMET_Phi_);
    
    // Pileup density
    outT->Branch("rho",&rho_);

    // Electron-positron vertex branches
    outT->Branch("nvtx",&nvtx_);
    outT->Branch("vtx_typ",&vtx_type_);
    outT->Branch("vtx_vxy", &vtx_recoVtxVxy_);
    outT->Branch("vtx_sigmavxy", &vtx_recoVtxSigmaVxy_);
    outT->Branch("vtx_vx",&vtx_recoVtxVx_);
    outT->Branch("vtx_vy",&vtx_recoVtxVy_);
    outT->Branch("vtx_vz",&vtx_recoVtxVz_);
    outT->Branch("vtx_reduced_chi2", &vtx_recoVtxReducedChi2_);
    outT->Branch("vtx_prob", &vtx_prob_);
    outT->Branch("vtx_dR",  &vtx_recoVtxDr_);
    outT->Branch("vtx_sign",&vtx_recoVtxSign_);
    outT->Branch("vtx_minDxy",&vtx_minDxy_);
    outT->Branch("vtx_METdPhi",&vtx_METdPhi_);
    outT->Branch("vtx_pt",&vtx_ll_pt_);
    outT->Branch("vtx_eta",&vtx_ll_eta_);
    outT->Branch("vtx_phi",&vtx_ll_phi_);
    outT->Branch("vtx_energy",&vtx_ll_e_);
    outT->Branch("vtx_m",&vtx_ll_m_);
    outT->Branch("vtx_px",&vtx_ll_px_);
    outT->Branch("vtx_py",&vtx_ll_py_);
    outT->Branch("vtx_pz",&vtx_ll_pz_);
    outT->Branch("vtx_PFIso4",&vtx_ll_PFIso_dR4_);
    outT->Branch("vtx_PFRelIso4",&vtx_ll_PFRelIso_dR4_);
    outT->Branch("vtx_PFIso3",&vtx_ll_PFIso_dR3_);
    outT->Branch("vtx_PFRelIso3",&vtx_ll_PFRelIso_dR3_);
    outT->Branch("vtx_PFIso8",&vtx_ll_PFIso_dR8_);
    outT->Branch("vtx_PFRelIso8",&vtx_ll_PFRelIso_dR8_);

    outT->Branch("vtx_e1_typ",&vtx_e1_type_);
    outT->Branch("vtx_e1_idx",&vtx_e1_idx_);
    outT->Branch("vtx_e2_typ",&vtx_e2_type_);
    outT->Branch("vtx_e2_idx",&vtx_e2_idx_);

    // Gen information
    if (!isData_) {
        outT->Branch("genWgt", &genwgt_);
        outT->Branch("genPU_obs", &genpuobs_);
        outT->Branch("genPU_true", &genputrue_);
        outT->Branch("nGenJet",&nGenJet_);
        outT->Branch("GenJet_pt", &genJetPt_);
        outT->Branch("GenJet_eta", &genJetEta_);
        outT->Branch("GenJet_phi", &genJetPhi_);
        outT->Branch("GenJet_METdPhi",&genJetMETdPhi_);
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
        outT->Branch("GenEle_vz",&genEleVz_);
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
        outT->Branch("GenPos_vz",&genPosVz_);
        outT->Branch("GenPosClosest_typ",&genPosClosestType_); // need to misspell type so it doesn't cause python issues
        outT->Branch("GenPosClosest_ind",&genPosClosestInd_);
        outT->Branch("GenPosClosest_dr",&genPosClosestDr_);
        outT->Branch("GenPosClosestReg_ind",&genPosClosestInd_reg_);
        outT->Branch("GenPosClosestReg_dr",&genPosClosestDr_reg_);
        outT->Branch("GenPosClosestLpt_ind",&genPosClosestInd_lpt_);
        outT->Branch("GenPosClosestLpt_dr",&genPosClosestDr_lpt_);

        // Gen Electron matches to isoTracks
        outT->Branch("nGenEleTrkMatches",&nGenEleTrkMatches);
        outT->Branch("genEleNearTk_pt",&genEleNearestTrack_pt);
        outT->Branch("genEleNearTk_eta",&genEleNearestTrack_eta);
        outT->Branch("genEleNearTk_phi",&genEleNearestTrack_phi);
        outT->Branch("genEleNearTk_dRgen",&genEleNearestTrack_dRGen);
        outT->Branch("genEleNearTk_pfIso3",&genEleNearestTrack_pfIso3);
        outT->Branch("genEleNearTk_miniIso",&genEleNearestTrack_miniIso);
        outT->Branch("genEleNearTk_dxy",&genEleNearestTrack_dxy);
        outT->Branch("genEleNearTk_dz",&genEleNearestTrack_dz);
        outT->Branch("genEleNearTk_highPurity",&genEleNearestTrack_highPurity);
        outT->Branch("genEleNearTk_loose",&genEleNearestTrack_Loose);
        outT->Branch("genEleNearTk_charge",&genEleNearestTrack_charge);
        outT->Branch("genEleNearTk_numPixHits",&genEleNearestTrack_numPixHits);
        outT->Branch("genEleNearTk_numStripHits",&genEleNearestTrack_numStripHits);
        outT->Branch("genEleNearTk_fromPV",&genEleNearestTrack_fromPV);
        outT->Branch("genEleNearTk_tkIdx",&genEleNearestTrack_tkIdx);

        // Gen Positron matches to isoTracks
        outT->Branch("nGenPosTrkMatches",&nGenPosTrkMatches);
        outT->Branch("genPosNearTk_pt",&genPosNearestTrack_pt);
        outT->Branch("genPosNearTk_eta",&genPosNearestTrack_eta);
        outT->Branch("genPosNearTk_phi",&genPosNearestTrack_phi);
        outT->Branch("genPosNearTk_dRgen",&genPosNearestTrack_dRGen);
        outT->Branch("genPosNearTk_pfIso3",&genPosNearestTrack_pfIso3);
        outT->Branch("genPosNearTk_miniIso",&genPosNearestTrack_miniIso);
        outT->Branch("genPosNearTk_dxy",&genPosNearestTrack_dxy);
        outT->Branch("genPosNearTk_dz",&genPosNearestTrack_dz);
        outT->Branch("genPosNearTk_highPurity",&genPosNearestTrack_highPurity);
        outT->Branch("genPosNearTk_loose",&genPosNearestTrack_Loose);
        outT->Branch("genPosNearTk_charge",&genPosNearestTrack_charge);
        outT->Branch("genPosNearTk_numPixHits",&genPosNearestTrack_numPixHits);
        outT->Branch("genPosNearTk_numStripHits",&genPosNearestTrack_numStripHits);
        outT->Branch("genPosNearTk_fromPV",&genPosNearestTrack_fromPV);
        outT->Branch("genPosNearTk_tkIdx",&genPosNearestTrack_tkIdx);

        // Gen Electron + Positron info
        outT->Branch("genEE_pt",&genEEPt_);
        outT->Branch("genEE_eta",&genEEEta_);
        outT->Branch("genEE_phi",&genEEPhi_);
        outT->Branch("genEE_energy",&genEEEn_);
        outT->Branch("genEE_mass",&genEEMass_);
        outT->Branch("genEE_dr",&genEEdR_);
        outT->Branch("genEE_METdPhi",&genEEMETdPhi_);
    }

}

void NtupleContainerV2::ClearTreeBranches() {
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
    genPosClosestType_ = -1;
    genPosClosestInd_ = -1;
    genPosClosestDr_ = 999.0;
    genPosClosestInd_reg_ = -1;
    genPosClosestDr_reg_ = 999.0;
    genPosClosestInd_lpt_ = -1;
    genPosClosestDr_lpt_ = 999.0;

    // Gen Electron matches to isoTracks
    nGenEleTrkMatches = 0;
    genEleNearestTrack_pt = -999.0;
    genEleNearestTrack_eta = -999.0;
    genEleNearestTrack_phi = -999.0;
    genEleNearestTrack_dRGen = -999.0;
    genEleNearestTrack_pfIso3 = -999.0;
    genEleNearestTrack_miniIso = -999.0;
    genEleNearestTrack_dxy = -999.0;
    genEleNearestTrack_dz = -999.0;
    genEleNearestTrack_highPurity = false;
    genEleNearestTrack_Loose = false;
    genEleNearestTrack_charge = 0;
    genEleNearestTrack_numPixHits = 0;
    genEleNearestTrack_numStripHits = 0;
    genEleNearestTrack_fromPV = 0;
    genEleNearestTrack_tkIdx = -1;

    // Gen Positron matches to isoTracks
    nGenPosTrkMatches = 0;
    genPosNearestTrack_pt = -999.0;
    genPosNearestTrack_eta = -999.0;
    genPosNearestTrack_phi = -999.0;
    genPosNearestTrack_dRGen = -999.0;
    genPosNearestTrack_pfIso3 = -999.0;
    genPosNearestTrack_miniIso = -999.0;
    genPosNearestTrack_dxy = -999.0;
    genPosNearestTrack_dz = -999.0;
    genPosNearestTrack_highPurity = false;
    genPosNearestTrack_Loose = false;
    genPosNearestTrack_charge = 0;
    genPosNearestTrack_numPixHits = 0;
    genPosNearestTrack_numStripHits = 0;
    genPosNearestTrack_fromPV = 0;
    genPosNearestTrack_tkIdx = -1;

    // Gen Electron + Positron info
    genEEPt_ = -999;
    genEEEta_ = -999;
    genEEPhi_ = -999;
    genEEEn_ = -999;
    genEEMass_ = -999;
    genEEdR_ = -999;
    genEEMETdPhi_ = -999;
    
    // Gen jet
    nGenJet_ = 0;
    genJetPt_.clear();
    genJetEta_.clear();
    genJetPhi_.clear();
    genJetMETdPhi_.clear();

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
    recoElectronID_cutVeto_.clear();
    recoElectronID_cutLoose_.clear();
    recoElectronID_cutMed_.clear();
    recoElectronID_cutTight_.clear();
    recoElectronID_mvaIso90_.clear();
    recoElectronID_mvaIso80_.clear();
    recoElectronID_mvaIsoLoose_.clear();
    recoElectronID_mva90_.clear();
    recoElectronID_mva80_.clear();
    recoElectronID_mvaLoose_.clear();
    recoElectronAngularRes_.clear();
    recoElectronE_.clear();
    recoElectronVxy_.clear();
    recoElectronVz_.clear();
    recoElectronDxy_.clear();
    recoElectronDxyError_.clear();
    recoElectronDz_.clear();
    recoElectronDzError_.clear();
    recoElectronTrkChi2_.clear();
    recoElectronTrkIso_.clear();
    recoElectronTrkRelIso_.clear();
    recoElectronCaloIso_.clear();
    recoElectronCaloRelIso_.clear();
    recoElectronPFIso_dR4_.clear();
    recoElectronPFRelIso_dR4_.clear();
    recoElectronPFIso_dR3_.clear();
    recoElectronPFRelIso_dR3_.clear();
    recoElectronPFIso_dR8_.clear();
    recoElectronPFRelIso_dR8_.clear();
    recoElectronPFIso_.clear();
    recoElectronPFRelIso_.clear();
    recoElectronChadIso_.clear();
    recoElectronNhadIso_.clear();
    recoElectronPhoIso_.clear();
    recoElectronRhoEA_.clear();
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
    recoLowPtElectronVxy_.clear();
    recoLowPtElectronVz_.clear();
    recoLowPtElectronDxy_.clear();
    recoLowPtElectronDxyError_.clear();
    recoLowPtElectronDz_.clear();
    recoLowPtElectronDzError_.clear();
    recoLowPtElectronTrkChi2_.clear();
    recoLowPtElectronTrkIso_.clear();
    recoLowPtElectronTrkRelIso_.clear();
    recoLowPtElectronCaloIso_.clear();
    recoLowPtElectronCaloRelIso_.clear();
    recoLowPtElectronPFIso_dR4_.clear();
    recoLowPtElectronPFRelIso_dR4_.clear();
    recoLowPtElectronPFIso_dR3_.clear();
    recoLowPtElectronPFRelIso_dR3_.clear();
    recoLowPtElectronPFIso_dR8_.clear();
    recoLowPtElectronPFRelIso_dR8_.clear();
    recoLowPtElectronPFIso_.clear();
    recoLowPtElectronPFRelIso_.clear();
    recoLowPtElectronChadIso_.clear();
    recoLowPtElectronNhadIso_.clear();
    recoLowPtElectronPhoIso_.clear();
    recoLowPtElectronRhoEA_.clear();
    recoLowPtElectronTrkProb_.clear();
    recoLowPtElectronTrkNumTrackerHits_.clear();
    recoLowPtElectronTrkNumPixHits_.clear();
    recoLowPtElectronTrkNumStripHits_.clear();
    recoLowPtElectronCharge_.clear();
    recoLowPtElectronMinDrToReg_.clear();

    // Gen weight and pileup
    genwgt_ = 0;
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
    PFNJetAll_ = 0;
    PFJetPt_.clear();
    PFJetEta_.clear();
    PFJetPhi_.clear();
    PFJetBTag_.clear();
    PFJetMETdPhi_.clear();
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
    PFMET_Pt_ = -999;
    PFMET_Phi_ = -999;
    PFMETJESUpPt_ = -9999;
    PFMETJESUpPhi_ = -9999;
    PFMETJESDownPt_ = -9999;
    PFMETJESDownPhi_ = -9999;
    PFMETJERUpPt_ = -9999;
    PFMETJERUpPhi_ = -9999;
    PFMETJERDownPt_ = -9999;
    PFMETJERDownPhi_ = -9999;
    PFMETJetResUpSmearPt_ = -9999;
    PFMETJetResUpSmearPhi_ = -9999;
    PFMETJetResDownSmearPt_ = -9999;
    PFMETJetResDownSmearPhi_ = -9999;
    
    CaloMET_ET_ = -999;
    CaloMET_Pt_ = -999;
    CaloMET_Phi_ = -999;

    // Pileup density
    rho_ = -9999;

    // Electron-positron vertices
    nvtx_ = 0;
    vtx_type_.clear();
    vtx_recoVtxVxy_.clear();
    vtx_recoVtxVz_.clear();
    vtx_recoVtxSigmaVxy_.clear();
    vtx_recoVtxVx_.clear();
    vtx_recoVtxVy_.clear();
    vtx_recoVtxVz_.clear();
    vtx_recoVtxReducedChi2_.clear();
    vtx_prob_.clear();
    vtx_recoVtxDr_.clear();
    vtx_recoVtxSign_.clear();
    vtx_minDxy_.clear();
    vtx_METdPhi_.clear();
    vtx_ll_pt_.clear();
    vtx_ll_eta_.clear();
    vtx_ll_phi_.clear();
    vtx_ll_e_.clear();
    vtx_ll_m_.clear();
    vtx_ll_px_.clear();
    vtx_ll_py_.clear();
    vtx_ll_pz_.clear();
    vtx_ll_PFIso_dR4_.clear();
    vtx_ll_PFRelIso_dR4_.clear();
    vtx_ll_PFIso_dR3_.clear();
    vtx_ll_PFRelIso_dR3_.clear();
    vtx_ll_PFIso_dR8_.clear();
    vtx_ll_PFRelIso_dR8_.clear();

    vtx_e1_type_.clear();
    vtx_e1_idx_.clear();
    vtx_e2_type_.clear();
    vtx_e2_idx_.clear();
}
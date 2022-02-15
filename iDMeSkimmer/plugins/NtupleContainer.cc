#include "NtupleContainer.hh"

NtupleContainer::NtupleContainer() : isData_(true) {}

NtupleContainer::~NtupleContainer() {}

void NtupleContainer::SetTree(TTree *tree, bool isData) { outT = tree; isData_ = isData; }

void NtupleContainer::CreateTreeBranches() {

    // Data or MC
    outT->Branch("isData",&isData_);

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
    outT->Branch("nEleVertexRegReg",&nEleVertex_regreg_);
    outT->Branch("EleVertexRegReg_vxy", &regreg_recoVtxVxy_);
    outT->Branch("EleVertexRegReg_vz",  &regreg_recoVtxVz_);
    outT->Branch("EleVertexRegReg_sigmavxy", &regreg_recoVtxSigmaVxy_);
    outT->Branch("EleVertexRegReg_x",&regreg_recoVtx_x_);
    outT->Branch("EleVertexRegReg_y",&regreg_recoVtx_y_);
    outT->Branch("EleVertexRegReg_z",&regreg_recoVtx_z_);
    outT->Branch("EleVertexRegReg_reduced_chi2", &regreg_recoVtxReducedChi2_);
    outT->Branch("EleVertexRegReg_dR",  &regreg_recoVtxDr_);
    outT->Branch("EleVertexRegReg_sign",&regreg_recoVtxSign_);
    outT->Branch("EleVertexRegReg_pt",&regreg_ll_pt_);
    outT->Branch("EleVertexRegReg_eta",&regreg_ll_eta_);
    outT->Branch("EleVertexRegReg_phi",&regreg_ll_phi_);
    outT->Branch("EleVertexRegReg_energy",&regreg_ll_e_);
    outT->Branch("EleVertexRegReg_m",&regreg_ll_m_);
    outT->Branch("EleVertexRegReg_px",&regreg_ll_px_);
    outT->Branch("EleVertexRegReg_py",&regreg_ll_py_);
    outT->Branch("EleVertexRegReg_pz",&regreg_ll_pz_);

    outT->Branch("nEleVertexLowLow",&nEleVertex_lowlow_);
    outT->Branch("EleVertexLowLow_vxy", &lowlow_recoVtxVxy_);
    outT->Branch("EleVertexLowLow_vz",  &lowlow_recoVtxVz_);
    outT->Branch("EleVertexLowLow_sigmavxy", &lowlow_recoVtxSigmaVxy_);
    outT->Branch("EleVertexLowLow_x",&lowlow_recoVtx_x_);
    outT->Branch("EleVertexLowLow_y",&lowlow_recoVtx_y_);
    outT->Branch("EleVertexLowLow_z",&lowlow_recoVtx_z_);
    outT->Branch("EleVertexLowLow_reduced_chi2", &lowlow_recoVtxReducedChi2_);
    outT->Branch("EleVertexLowLow_dR",  &lowlow_recoVtxDr_);
    outT->Branch("EleVertexLowLow_sign",&lowlow_recoVtxSign_);
    outT->Branch("EleVertexLowLow_pt",&lowlow_ll_pt_);
    outT->Branch("EleVertexLowLow_eta",&lowlow_ll_eta_);
    outT->Branch("EleVertexLowLow_phi",&lowlow_ll_phi_);
    outT->Branch("EleVertexLowLow_energy",&lowlow_ll_e_);
    outT->Branch("EleVertexLowLow_m",&lowlow_ll_m_);
    outT->Branch("EleVertexLowLow_px",&lowlow_ll_px_);
    outT->Branch("EleVertexLowLow_py",&lowlow_ll_py_);
    outT->Branch("EleVertexLowLow_pz",&lowlow_ll_pz_);

    outT->Branch("nEleVertexLowReg",&nEleVertex_lowreg_);
    outT->Branch("EleVertexLowReg_vxy", &lowreg_recoVtxVxy_);
    outT->Branch("EleVertexLowReg_vz",  &lowreg_recoVtxVz_);
    outT->Branch("EleVertexLowReg_sigmavxy", &lowreg_recoVtxSigmaVxy_);
    outT->Branch("EleVertexLowReg_x",&lowreg_recoVtx_x_);
    outT->Branch("EleVertexLowReg_y",&lowreg_recoVtx_y_);
    outT->Branch("EleVertexLowReg_z",&lowreg_recoVtx_z_);
    outT->Branch("EleVertexLowReg_reduced_chi2", &lowreg_recoVtxReducedChi2_);
    outT->Branch("EleVertexLowReg_dR",  &lowreg_recoVtxDr_);
    outT->Branch("EleVertexLowReg_sign",&lowreg_recoVtxSign_);
    outT->Branch("EleVertexLowReg_pt",&lowreg_ll_pt_);
    outT->Branch("EleVertexLowReg_eta",&lowreg_ll_eta_);
    outT->Branch("EleVertexLowReg_phi",&lowreg_ll_phi_);
    outT->Branch("EleVertexLowReg_energy",&lowreg_ll_e_);
    outT->Branch("EleVertexLowReg_m",&lowreg_ll_m_);
    outT->Branch("EleVertexLowReg_px",&lowreg_ll_px_);
    outT->Branch("EleVertexLowReg_py",&lowreg_ll_py_);
    outT->Branch("EleVertexLowReg_pz",&lowreg_ll_pz_);

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
        outT->Branch("GenPart_vz", &genVz_);
        outT->Branch("GenPart_x",&genVtx_x_);
        outT->Branch("GenPart_y",&genVtx_y_);
        outT->Branch("GenPart_z",&genVtx_z_);
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
    genVz_.clear();
    genVtx_x_.clear();
    genVtx_y_.clear();
    genVtx_z_.clear();
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
    nEleVertex_regreg_ = 0;
    regreg_recoVtxVxy_.clear();
    regreg_recoVtxVz_.clear();
    regreg_recoVtxSigmaVxy_.clear();
    regreg_recoVtx_x_.clear();
    regreg_recoVtx_y_.clear();
    regreg_recoVtx_z_.clear();
    regreg_recoVtxReducedChi2_.clear();
    regreg_recoVtxDr_.clear();
    regreg_recoVtxSign_.clear();
    regreg_ll_pt_.clear();
    regreg_ll_eta_.clear();
    regreg_ll_phi_.clear();
    regreg_ll_e_.clear();
    regreg_ll_m_.clear();
    regreg_ll_px_.clear();
    regreg_ll_py_.clear();
    regreg_ll_pz_.clear();

    nEleVertex_lowlow_ = 0;
    lowlow_recoVtxVxy_.clear();
    lowlow_recoVtxVz_.clear();
    lowlow_recoVtxSigmaVxy_.clear();
    lowlow_recoVtx_x_.clear();
    lowlow_recoVtx_y_.clear();
    lowlow_recoVtx_z_.clear();
    lowlow_recoVtxReducedChi2_.clear();
    lowlow_recoVtxDr_.clear();
    lowlow_recoVtxSign_.clear();
    lowlow_ll_pt_.clear();
    lowlow_ll_eta_.clear();
    lowlow_ll_phi_.clear();
    lowlow_ll_e_.clear();
    lowlow_ll_m_.clear();
    lowlow_ll_px_.clear();
    lowlow_ll_py_.clear();
    lowlow_ll_pz_.clear();

    nEleVertex_lowreg_ = 0;
    lowreg_recoVtxVxy_.clear();
    lowreg_recoVtxVz_.clear();
    lowreg_recoVtxSigmaVxy_.clear();
    lowreg_recoVtx_x_.clear();
    lowreg_recoVtx_y_.clear();
    lowreg_recoVtx_z_.clear();
    lowreg_recoVtxReducedChi2_.clear();
    lowreg_recoVtxDr_.clear();
    lowreg_recoVtxSign_.clear();
    lowreg_ll_pt_.clear();
    lowreg_ll_eta_.clear();
    lowreg_ll_phi_.clear();
    lowreg_ll_e_.clear();
    lowreg_ll_m_.clear();
    lowreg_ll_px_.clear();
    lowreg_ll_py_.clear();
    lowreg_ll_pz_.clear();
}
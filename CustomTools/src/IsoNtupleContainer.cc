#include "iDMe/CustomTools/interface/IsoNtupleContainer.hh"

IsoNtupleContainer::IsoNtupleContainer() {}

IsoNtupleContainer::~IsoNtupleContainer() {}

void IsoNtupleContainer::CreateTreeBranches() {

    // Create branches from parent class
    NtupleContainer::CreateTreeBranches();

    // PF Candidates
    outT->Branch("nPF",&nPFCand_);
    outT->Branch("PF_ID",&pfID_);
    outT->Branch("PF_charge",&pfCharge_);
    outT->Branch("PF_pt",&pfPt_);
    outT->Branch("PF_eta",&pfEta_);
    outT->Branch("PF_phi",&pfPhi_);
    outT->Branch("PF_E",&pfE_);

}

void IsoNtupleContainer::ClearTreeBranches() {

    // Calling parent function
    NtupleContainer::ClearTreeBranches();

    // PF Candidates
    nPFCand_ = 0;
    pfID_.clear();
    pfCharge_.clear();
    pfPt_.clear();
    pfEta_.clear();
    pfPhi_.clear();
    pfE_.clear();
    
}
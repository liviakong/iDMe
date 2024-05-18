#include "iDMe/CustomTools/interface/IsolationCalculator.hh"

IsolationCalculator::~IsolationCalculator() {}

IsolationCalculator::IsolationCalculator(edm::Handle<vector<pat::Electron> > &recoElectronHandle_, edm::Handle<vector<pat::Electron> > &lowPtElectronHandle_, edm::Handle<vector<pat::PackedCandidate> > PFcandHandle_, NtupleContainerV2 &nt_)
:
electronHandle(recoElectronHandle_), lowPtElectronHandle(lowPtElectronHandle_), PFcandHandle(PFcandHandle_), nt(nt_)
{}

void IsolationCalculator::calcIso() {
    /*nt.recoElectronPFIso_dR4_.resize((*electronHandle).size(),0.0);
    nt.recoElectronPFRelIso_dR4_.resize((*electronHandle).size(),0.0);
    nt.recoElectronPFIso_dR3_.resize((*electronHandle).size(),0.0);
    nt.recoElectronPFRelIso_dR3_.resize((*electronHandle).size(),0.0);
    nt.recoElectronPFIso_dR8_.resize((*electronHandle).size(),0.0);
    nt.recoElectronPFRelIso_dR8_.resize((*electronHandle).size(),0.0);

    nt.recoLowPtElectronPFIso_dR4_.resize((*lowPtElectronHandle).size(),0.0);
    nt.recoLowPtElectronPFRelIso_dR4_.resize((*lowPtElectronHandle).size(),0.0);
    nt.recoLowPtElectronPFIso_dR3_.resize((*lowPtElectronHandle).size(),0.0);
    nt.recoLowPtElectronPFRelIso_dR3_.resize((*lowPtElectronHandle).size(),0.0);
    nt.recoLowPtElectronPFIso_dR8_.resize((*lowPtElectronHandle).size(),0.0);
    nt.recoLowPtElectronPFRelIso_dR8_.resize((*lowPtElectronHandle).size(),0.0);

    nt.vtx_ll_PFIso_dR4_.resize(nt.nvtx_,0.0);
    nt.vtx_ll_PFRelIso_dR4_.resize(nt.nvtx_,0.0);
    nt.vtx_ll_PFIso_dR3_.resize(nt.nvtx_,0.0);
    nt.vtx_ll_PFRelIso_dR3_.resize(nt.nvtx_,0.0);
    nt.vtx_ll_PFIso_dR8_.resize(nt.nvtx_,0.0);
    nt.vtx_ll_PFRelIso_dR8_.resize(nt.nvtx_,0.0);

    for (auto & pfc : *PFcandHandle) {
        // Calculate electron isolation
        for (unsigned int i = 0; i < (*electronHandle).size(); i++) {
            auto ele = (*electronHandle)[i];
            float dR = reco::deltaR(ele,pfc);
            if (dR < 0.015) continue; // veto if too close
            float pt = pfc.pt();
            float rel = pfc.pt()/ele.pt();
            if (dR < 0.8) {
                nt.recoElectronPFIso_dR8_[i] += pt;
                nt.recoElectronPFRelIso_dR8_[i] += rel;
            }
            if (dR < 0.4) {
                nt.recoElectronPFIso_dR4_[i] += pt;
                nt.recoElectronPFRelIso_dR4_[i] += rel;
            }
            if (dR < 0.3) {
                nt.recoElectronPFIso_dR3_[i] += pt;
                nt.recoElectronPFRelIso_dR3_[i] += rel;
            }
        }

        // Calculate low pT electron isolation
        for (unsigned int i = 0; i < (*lowPtElectronHandle).size(); i++) {
            auto ele = (*lowPtElectronHandle)[i];
            float dR = reco::deltaR(ele,pfc);
            if (dR < 0.015) continue; // veto if too close
            float pt = pfc.pt();
            float rel = pfc.pt()/ele.pt();
            if (dR < 0.8) {
                nt.recoLowPtElectronPFIso_dR8_[i] += pt;
                nt.recoLowPtElectronPFRelIso_dR8_[i] += rel;
            }
            if (dR < 0.4) {
                nt.recoLowPtElectronPFIso_dR4_[i] += pt;
                nt.recoLowPtElectronPFRelIso_dR4_[i] += rel;
            }
            if (dR < 0.3) {
                nt.recoLowPtElectronPFIso_dR3_[i] += pt;
                nt.recoLowPtElectronPFRelIso_dR3_[i] += rel;
            }
        }

        // Calculate vertex isolation
        for (int i = 0; i < nt.nvtx_; i++) {
            int i1 = nt.vtx_e1_idx_[i];
            auto e1 = nt.vtx_e1_type_[i] == "R" ? (*electronHandle)[i1] : (*lowPtElectronHandle)[i1];
            int i2 = nt.vtx_e2_idx_[i];
            auto e2 = nt.vtx_e2_type_[i] == "R" ? (*electronHandle)[i2] : (*lowPtElectronHandle)[i2];
            if ((reco::deltaR(e1,pfc) < 0.015) || (reco::deltaR(e2,pfc) < 0.015)) continue;


            auto ll = e1.p4() + e2.p4();
            float pt = pfc.pt();
            float rel = pfc.pt()/ll.pt();
            float dR = reco::deltaR(ll,pfc);
            if (dR < 0.8) {
                nt.vtx_ll_PFIso_dR8_[i] += pt;
                nt.vtx_ll_PFRelIso_dR8_[i] += rel;
            }
            if (dR < 0.4) {
                nt.vtx_ll_PFIso_dR4_[i] += pt;
                nt.vtx_ll_PFRelIso_dR4_[i] += rel;
            }
            if (dR < 0.3) {
                nt.vtx_ll_PFIso_dR3_[i] += pt;
                nt.vtx_ll_PFRelIso_dR3_[i] += rel;
            }
        }
    }*/
}
#include "iDMe/CustomTools/interface/Helpers.hh"

Helper::Helper() {}
Helper::~Helper() {}

bool Helper::JetID(const pat::Jet &jet, int year) {
    auto eta = jet.eta();
    auto neutEmFrac = jet.neutralEmEnergyFraction();
    auto neutHadFrac = jet.neutralHadronEnergyFraction();
    auto nConstit = jet.nConstituents();
    auto chargedHadFrac = jet.chargedHadronEnergyFraction();
    auto chargedMult = jet.chargedMultiplicity();
    auto neutMult = jet.neutralMultiplicity();

    bool passID = false;

    if (year == 2016) {
        if (abs(eta) <= 2.4) {
            passID = (neutHadFrac < 0.9) && (neutEmFrac < 0.9) && (nConstit > 1) && (chargedHadFrac > 0) && (chargedMult > 0);
        }
        else if ((abs(eta) >= 2.4) && (abs(eta) <= 2.7)) {
            passID = (neutHadFrac < 0.9) && (neutEmFrac < 0.99);
        }
        else if ((abs(eta) >= 2.7) && (abs(eta) <= 3.0)) {
            passID = (neutHadFrac < 0.9) && (neutEmFrac < 0.99) && (neutEmFrac > 0) && (neutMult > 1);
        }
        else if ((abs(eta) >= 3.0) && (abs(eta) < 5.0)) {
            passID = (neutHadFrac > 0.2) && (neutEmFrac < 0.9) && (neutMult > 10);
        }
    }
    else if ((year == 2017) || (year == 2018)) {
        if (abs(eta) <= 2.6) {
            passID = (neutHadFrac < 0.9) && (neutEmFrac < 0.9) && (nConstit > 1) && (chargedHadFrac > 0) && (chargedMult > 0);
        }
        else if ((abs(eta) >= 2.6) && (abs(eta) <= 2.7)) {
            passID = (neutHadFrac < 0.9) && (neutEmFrac < 0.99) && (chargedMult > 0);
        }
        else if ((abs(eta) >= 2.7) && (abs(eta) <= 3.0)) {
            passID = (neutEmFrac < 0.99) && (neutEmFrac > 0.01) && (neutMult > 1);
        }
        else if ((abs(eta) >= 3.0) && (abs(eta) < 5.0)) {
            passID = (neutHadFrac > 0.2) && (neutEmFrac < 0.9) && (neutMult > 10);
        }
    }

    // Apply additional cuts to leading jet (see monojet analysis: https://arxiv.org/pdf/1703.01651.pdf)
    //if (idx == 0) {
    //    passID = passID && (chargedHadFrac > 0.1) && (neutHadFrac < 0.8);
    //}

    return passID;
}
#ifndef JETCORRECTIONS_HH
#define JETCORRECTIONS_HH

#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>
#include <boost/any.hpp>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "iDMeAnalysis/CustomTools/interface/NtupleContainer.hh"

using std::vector;

class JetCorrections {
    public:
        // methods
        JetCorrections(edm::Handle<vector<reco::PFJet> > &recoJetHandle_, edm::Handle<reco::JetCorrector> &jetCorrectorHandle_, NtupleContainer &nt_, edm::Handle<reco::JetTagCollection> &bTagProbbHandle_, edm::Handle<reco::JetTagCollection> &bTagProbbbHandle_, edm::Handle<vector<reco::GsfElectron> > recoElectronHandle_, int year_, bool isData_);
        
        void Correct(JME::JetResolution &resolution, JME::JetResolutionScaleFactor &resolution_sf, JetCorrectionUncertainty &jecUnc, edm::Handle<double> &rhoHandle, edm::Handle<vector<reco::GenJet> > &genJetHandle, reco::PFMET &PFMET);
        
        bool passJetID(const reco::PFJet &jet, int idx);
        
        ~JetCorrections();

        // variables
        edm::Handle<vector<reco::PFJet> > &jetHandle;
        edm::Handle<reco::JetCorrector> &corrHandle;
        edm::Handle<vector<reco::GsfElectron> > &eleHandle;
        NtupleContainer &nt;
        edm::Handle<reco::JetTagCollection> &bTagProbbHandle;
        edm::Handle<reco::JetTagCollection> &bTagProbbbHandle;
        int year;
        bool isData;
        std::mt19937 m_random_generator;

        float corr_METpx, corr_METpy;
        float corr_METpx_JESUp, corr_METpy_JESUp;
        float corr_METpx_JESDown, corr_METpy_JESDown;
        float corr_METpx_JERUp, corr_METpy_JERUp;
        float corr_METpx_JERDown, corr_METpy_JERDown;
};

#endif
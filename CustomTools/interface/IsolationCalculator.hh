#ifndef ISOLATIONCALCULATOR_HH
#define ISOLATIONCALCULATOR_HH

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

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "iDMe/CustomTools/interface/NtupleContainerV2.hh"

using std::vector;

class IsolationCalculator {
    public:
        // methods
        IsolationCalculator(edm::Handle<vector<pat::Electron> > &recoElectronHandle_, edm::Handle<vector<pat::Electron> > &lowPtElectronHandle_, edm::Handle<vector<pat::PackedCandidate> > PFcandHandle_, NtupleContainerV2 &nt_);
        
        void calcIso();
        
        ~IsolationCalculator();

        // variables
        edm::Handle<vector<pat::Electron> > electronHandle;
        edm::Handle<vector<pat::Electron> > lowPtElectronHandle;
        edm::Handle<vector<pat::PackedCandidate> > PFcandHandle;
        NtupleContainerV2 &nt;
};

#endif
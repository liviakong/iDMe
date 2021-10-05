// -*- C++ -*-
//
// Package:    iDMeAnalysis/iDMeSkimmer
// Class:      iDMeSkimmer
//
/**\class iDMeSkimmer iDMeSkimmer.cc iDMeAnalysis/iDMeSkimmer/plugins/iDMeSkimmer.cc

 Description: MiniAOD skimmer for iDM analysis with electrons
*/
//
// Original Author:  Samuel Bright-Thonney
//         Created:  Tue, 21 Sep 2021 17:00:38 GMT
//

#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>
#include <boost/any.hpp>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"

#include "NtupleContainer.hh"

class iDMeSkimmer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit iDMeSkimmer(const edm::ParameterSet&);
      ~iDMeSkimmer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      bool getCollections(const edm::Event&);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      TTree *recoT, *genT;
      NtupleContainer nt;
      edm::Service<TFileService> fs;

      std::mt19937 m_random_generator;

      bool isData;

      // Tokens 
      const edm::EDGetTokenT<pat::ElectronCollection> recoElectronToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> lowPtElectronToken_;
      const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
      const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
      const edm::EDGetTokenT<reco::GenMETCollection> genMETToken_;
      const edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;

      // Handles
      edm::Handle<pat::ElectronCollection> recoElectronHandle_;
      edm::Handle<pat::ElectronCollection> lowPtElectronHandle_;
      edm::Handle<reco::GenParticleCollection> genParticleHandle_;
      edm::Handle<reco::GenJetCollection> genJetHandle_;
      edm::Handle<reco::GenMETCollection> genMETHandle_;
      edm::Handle<GenEventInfoProduct> genEvtInfoHandle_;
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
iDMeSkimmer::iDMeSkimmer(const edm::ParameterSet& ps)
 :
   isData(ps.getParameter<bool>("isData")),
   recoElectronToken_(consumes<pat::ElectronCollection>(ps.getParameter<edm::InputTag>("recoElectron"))),
   lowPtElectronToken_(consumes<pat::ElectronCollection>(ps.getParameter<edm::InputTag>("lowPtElectron"))),
   genParticleToken_(consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("genParticle"))),
   genJetToken_(consumes<reco::GenJetCollection>(ps.getParameter<edm::InputTag>("genJet"))),
   genMETToken_(consumes<reco::GenMETCollection>(ps.getParameter<edm::InputTag>("genMET"))),
   genEvtInfoToken_(consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("genEvt")))
{
   usesResource("TFileService");
   m_random_generator = std::mt19937(37428479);

}


iDMeSkimmer::~iDMeSkimmer() = default;


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void
iDMeSkimmer::beginJob()
{
   recoT = fs->make<TTree>("recoT", "recoT");
   nt.SetRecoTree(recoT);
   if (!isData) {
      genT = fs->make<TTree>("genT", "genT");
      nt.SetGenTree(genT);
   }
   nt.CreateTreeBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void
iDMeSkimmer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
iDMeSkimmer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   edm::ParameterSetDescription desc;

   desc.add<bool>("isData", 0);

   desc.add<edm::InputTag>("recoElectron",edm::InputTag("slimmedElectrons"));
   desc.add<edm::InputTag>("lowPtElectron",edm::InputTag("slimmedLowPtElectrons"));
   desc.add<edm::InputTag>("genParticle",edm::InputTag("prunedGenParticles"));
   desc.add<edm::InputTag>("genJet",edm::InputTag("slimmedGenJets"));
   desc.add<edm::InputTag>("genMET",edm::InputTag("genMetTrue"));
   desc.add<edm::InputTag>("genEvt", edm::InputTag("generator"));
   
   descriptions.add("iDMeSkimmer", desc);
}

// ------------ method called for each event  ------------
void
iDMeSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using std::cout, std::endl; 

   // Retrieving event data and assigning to handles
   iEvent.getByToken(recoElectronToken_,recoElectronHandle_);
   iEvent.getByToken(lowPtElectronToken_,lowPtElectronHandle_);
   iEvent.getByToken(genParticleToken_,genParticleHandle_);
   iEvent.getByToken(genJetToken_,genJetHandle_);
   iEvent.getByToken(genMETToken_,genMETHandle_);
   iEvent.getByToken(genEvtInfoToken_,genEvtInfoHandle_);

   // Clear tree branches before filling
   nt.ClearTreeBranches();

   //////
   // Computing derived quantities and filling the trees
   //////
   nt.eventNum_ = iEvent.id().event();
   nt.lumiSec_ = iEvent.luminosityBlock();
   nt.runNum_ = iEvent.id().run();
   
   // Handling default electrons
   nt.nElectronDefault_ = recoElectronHandle_->size();
   std::vector<reco::TrackRef> reg_eleTracks{};
   for (unsigned int i = 0; i < recoElectronHandle_->size(); i++) {
      pat::ElectronRef ele(recoElectronHandle_,i);
      nt.recoElectronPt_.push_back(ele->pt());
      nt.recoElectronEta_.push_back(ele->eta());
      nt.recoElectronPhi_.push_back(ele->phi());
      nt.recoElectronVxy_.push_back(ele->trackPositionAtVtx().rho());
      nt.recoElectronVz_.push_back(ele->trackPositionAtVtx().z());
      nt.recoElectronCharge_.push_back(ele->charge());
      // Filling tracks
      auto track = ele->closestCtfTrackRef();
         if (!track.isNonnull())
            cout << "Track " << i << " from regular ele reco is not valid! " << endl;
         else {
            reg_eleTracks.emplace_back(track);
         }
   }

   // Handling low-pT electrons
   nt.nElectronLowPt_ = lowPtElectronHandle_->size();
   std::vector<reco::TrackRef> lowpt_eleTracks{};
   for (unsigned int i = 0; i < lowPtElectronHandle_->size(); i++) {
      pat::ElectronRef ele(lowPtElectronHandle_,i);
      nt.recoLowPtElectronPt_.push_back(ele->pt());
      nt.recoLowPtElectronPhi_.push_back(ele->phi());
      nt.recoLowPtElectronEta_.push_back(ele->eta());
      nt.recoLowPtElectronVxy_.push_back(ele->trackPositionAtVtx().rho());
      nt.recoLowPtElectronVz_.push_back(ele->trackPositionAtVtx().z());
      nt.recoLowPtElectronCharge_.push_back(ele->charge());
      // Filling tracks
      auto track = ele->closestCtfTrackRef();
         if (!track.isNonnull())
            cout << "Track " << i << " from low-pT ele reco is not valid! " << endl;
         else {
            lowpt_eleTracks.emplace_back(track);
         }
   }

   // Reconstructing electron vertices
   // Sort electron tracks by pT
   std::sort(reg_eleTracks.begin(), reg_eleTracks.end(), [](const auto & l, const auto & r) {
            return l->pt() > r->pt();
            });
   std::sort(lowpt_eleTracks.begin(), lowpt_eleTracks.end(), [](const auto & l, const auto & r) {
            return l->pt() > r->pt();
            });
   // Set up objects for vertex reco
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   KalmanVertexFitter kvf(true);
   // Define vertex reco function 
   auto computeVertices = [&](vector<reco::TrackRef> coll_1, vector<reco::TrackRef> coll_2, std::string type) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
               reco::TrackRef ele_i, ele_j;
               if (i < coll_1.size()) {
                  ele_i = coll_1[i];
               }
               if (j < coll_2.size()) {
                  ele_j = coll_2[j];
               }

               TransientVertex tv;
               if (ele_i.isNonnull() && ele_j.isNonnull() && i != j) {
                  vector<reco::TransientTrack> transient_tracks{};
                  transient_tracks.push_back(theB->build(ele_i));
                  transient_tracks.push_back(theB->build(ele_j));
                  tv = kvf.vertex(transient_tracks);
               }
               float vxy = -9999;
               float sigma_vxy = -9999;
               float vtx_chi2 = 999999;
               float vz = -9999;
               float dr = -9999;
               if (tv.isValid()) {
                  reco::Vertex vertex = reco::Vertex(tv);
                  vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
                  sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                           vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
                  //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
                  vtx_chi2 = vertex.normalizedChi2();
                  vz = vertex.z();
                  dr = reco::deltaR(*ele_i, *ele_j);
               }

               if (type == "lowlow") {
                  nt.lowlow_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.lowlow_recoVtxVxy_.push_back(vxy);
                  nt.lowlow_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.lowlow_recoVtxVz_.push_back(vz);
                  nt.lowlow_recoVtxDr_.push_back(dr);
               }
               else if (type == "regreg") {
                  nt.regreg_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.regreg_recoVtxVxy_.push_back(vxy);
                  nt.regreg_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.regreg_recoVtxVz_.push_back(vz);
                  nt.regreg_recoVtxDr_.push_back(dr);
               }
               else if (type == "lowreg") {
                  nt.lowreg_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.lowreg_recoVtxVxy_.push_back(vxy);
                  nt.lowreg_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.lowreg_recoVtxVz_.push_back(vz);
                  nt.lowreg_recoVtxDr_.push_back(dr);
               }
            }
        }
    };

    // lowpT-lowpT
    computeVertices(lowpt_eleTracks, lowpt_eleTracks, "lowlow");
    // regular-regular
    computeVertices(reg_eleTracks, reg_eleTracks, "regreg");
    // lowpT-regular
    computeVertices(lowpt_eleTracks, reg_eleTracks, "lowreg");

   ///////////////////////////////////////

   //Handling gen particles
   
   if (!isData) {
      nt.nGen_ = (int)genParticleHandle_->size();
      // Gen weight
      nt.genwgt_ = genEvtInfoHandle_->weight();
      for (size_t i = 0; i < genParticleHandle_->size(); i++) {
         reco::GenParticleRef genParticle(genParticleHandle_, i);
         if (!genParticle->isHardProcess()) continue;
         nt.genID_.push_back(genParticle->pdgId());
         nt.genCharge_.push_back(genParticle->charge());
         nt.genPt_.push_back(genParticle->pt());
         nt.genEta_.push_back(genParticle->eta());
         nt.genPhi_.push_back(genParticle->phi());
         nt.genPz_.push_back(genParticle->pz());
         nt.genEn_.push_back(genParticle->energy());
         nt.genVxy_.push_back(genParticle->vertex().rho());
         nt.genVz_.push_back(genParticle->vz());
         nt.genMass_.push_back(genParticle->mass());
      }
      
      // all gen jets
      for (size_t i = 0; i < genJetHandle_->size(); i++) {
         reco::GenJetRef jetRef(genJetHandle_, i);
         nt.genJetPt_.push_back(jetRef->pt());
         nt.genJetEta_.push_back(jetRef->eta());
         nt.genJetPhi_.push_back(jetRef->phi());
      }

      // Lead gen MET
      if (genMETHandle_->size() > 0) {
         reco::GenMETRef metRef(genMETHandle_, 0);
         nt.genLeadMETPt_ = metRef->pt();
         nt.genLeadMETPhi_ = metRef->phi();
      }

      genT->Fill();
    }


   recoT->Fill();
   return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(iDMeSkimmer);

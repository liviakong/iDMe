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
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
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
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

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
      TTree *outT;
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
      const edm::EDGetTokenT<reco::VertexCollection> primaryVertexToken_;
      const edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
      const edm::EDGetTokenT<vector<pat::Photon> > photonsToken_;
      const edm::EDGetTokenT<vector<pat::MET> > METToken_;
      const edm::EDGetTokenT<vector<pat::MET> > PuppiMETToken_;

      // Handles
      edm::Handle<pat::ElectronCollection> recoElectronHandle_;
      edm::Handle<pat::ElectronCollection> lowPtElectronHandle_;
      edm::Handle<reco::GenParticleCollection> genParticleHandle_;
      edm::Handle<reco::GenJetCollection> genJetHandle_;
      edm::Handle<reco::GenMETCollection> genMETHandle_;
      edm::Handle<GenEventInfoProduct> genEvtInfoHandle_;
      edm::Handle<reco::VertexCollection> primaryVertexHandle_;
      edm::Handle<vector<reco::Conversion> > conversionsHandle_;
      edm::Handle<vector<pat::Photon> > photonsHandle_;
      edm::Handle<vector<pat::MET> > METHandle_;
      edm::Handle<vector<pat::MET> > PuppiMETHandle_;
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
   genEvtInfoToken_(consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("genEvt"))),
   primaryVertexToken_(consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("primaryVertex"))),
   conversionsToken_(consumes<vector<reco::Conversion> >(ps.getParameter<edm::InputTag>("conversions"))),
   photonsToken_(consumes<vector<pat::Photon> >(ps.getParameter<edm::InputTag>("photons"))),
   METToken_(consumes<vector<pat::MET> >(ps.getParameter<edm::InputTag>("MET"))),
   PuppiMETToken_(consumes<vector<pat::MET> >(ps.getParameter<edm::InputTag>("PuppiMET")))
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
   outT = fs->make<TTree>("outT", "outT");
   nt.SetTree(outT,isData);
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
   desc.add<edm::InputTag>("primaryVertex",edm::InputTag("offlineSlimmedPrimaryVertices"));
   desc.add<edm::InputTag>("conversions",edm::InputTag("reducedEgamma","reducedConversions","PAT"));
   desc.add<edm::InputTag>("photons",edm::InputTag("slimmedPhotons"));
   desc.add<edm::InputTag>("MET",edm::InputTag("slimmedMETs"));
   desc.add<edm::InputTag>("PuppiMET",edm::InputTag("slimmedMETsPuppi"));
   
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
   iEvent.getByToken(primaryVertexToken_,primaryVertexHandle_);
   iEvent.getByToken(conversionsToken_,conversionsHandle_);
   iEvent.getByToken(conversionsToken_,conversionsHandle_);
   iEvent.getByToken(METToken_,METHandle_);
   iEvent.getByToken(PuppiMETToken_,PuppiMETHandle_);

   // Clear tree branches before filling
   nt.ClearTreeBranches();

   //////
   // Computing derived quantities and filling the trees
   //////
   nt.eventNum_ = iEvent.id().event();
   nt.lumiSec_ = iEvent.luminosityBlock();
   nt.runNum_ = iEvent.id().run();
   reco::Vertex pv = *primaryVertexHandle_->begin();

   // Handling Regular MET
   bool foundPF = false; bool foundCalo = false;
   for (const auto & met : *METHandle_) {
      // Take the first PF and Calo MET entries to be the ones we want to use (there shouldn't be > 1 of each?)
      // Recording "type 1" MET with the cor[VAR] functions (default in miniAOD)
      if (foundPF && foundCalo) break;
      if (met.isPFMET() && !foundPF) {
         foundPF = true;
         nt.PFMET_ET_ = met.corSumEt();
         nt.PFMET_Px_ = met.corPx();
         nt.PFMET_Py_ = met.corPy();
         nt.PFMET_Pt_ = met.corPt();
         nt.PFMET_Phi_ = met.corPhi();
      }
      if (met.isCaloMET() && !foundCalo) {
         foundCalo = true;
         nt.CaloMET_ET_ = met.corSumEt();
         nt.CaloMET_Px_ = met.corPx();
         nt.CaloMET_Py_ = met.corPy();
         nt.CaloMET_Pt_ = met.corPt();
         nt.CaloMET_Phi_ = met.corPhi();
      }
   }

   // Handling Puppi MET
   foundPF = false; foundCalo = false;
   for (const auto & met : *PuppiMETHandle_) {
      // Take the first PF and Calo MET entries to be the ones we want to use (there shouldn't be > 1 of each?)
      // Recording "type 1" MET with the cor[VAR] functions (default in miniAOD)
      if (foundPF && foundCalo) break;
      if (met.isPFMET() && !foundPF) {
         foundPF = true;
         nt.PuppiPFMET_ET_ = met.corSumEt();
         nt.PuppiPFMET_Px_ = met.corPx();
         nt.PuppiPFMET_Py_ = met.corPy();
         nt.PuppiPFMET_Pt_ = met.corPt();
         nt.PuppiPFMET_Phi_ = met.corPhi();
      }
      if (met.isCaloMET() && !foundCalo) {
         foundCalo = true;
         nt.PuppiCaloMET_ET_ = met.corSumEt();
         nt.PuppiCaloMET_Px_ = met.corPx();
         nt.PuppiCaloMET_Py_ = met.corPy();
         nt.PuppiCaloMET_Pt_ = met.corPt();
         nt.PuppiCaloMET_Phi_ = met.corPhi();
      }
   }
   
   // Handling default electrons
   nt.nElectronDefault_ = recoElectronHandle_->size();
   std::vector<reco::GsfTrackRef> reg_eleTracks{};
   for (unsigned int i = 0; i < recoElectronHandle_->size(); i++) {
      pat::ElectronRef ele(recoElectronHandle_,i);
      nt.recoElectronPt_.push_back(ele->pt());
      nt.recoElectronEta_.push_back(ele->eta());
      nt.recoElectronPhi_.push_back(ele->phi());
      nt.recoElectronE_.push_back(ele->energy());
      nt.recoElectronPx_.push_back(ele->px());
      nt.recoElectronPy_.push_back(ele->py());
      nt.recoElectronPz_.push_back(ele->pz());
      nt.recoElectronVxy_.push_back(ele->trackPositionAtVtx().rho());
      nt.recoElectronVz_.push_back(ele->trackPositionAtVtx().z());
      nt.recoElectronTrkIso_.push_back(ele->trackIso());
      nt.recoElectronCharge_.push_back(ele->charge());
      nt.recoElectron_passConversionVeto_.push_back((int)ele->passConversionVeto());
      // Filling tracks
      reco::GsfTrackRef track = ele->gsfTrack();
      if (!track.isNonnull())
         cout << "Track " << i << " from regular ele reco is not valid! " << endl;
      else {
         reg_eleTracks.emplace_back(track);
         nt.recoElectronDxy_.push_back(track->dxy(pv.position()));
         nt.recoElectronDxyError_.push_back(track->dxyError());
         nt.recoElectronDz_.push_back(track->dz(pv.position()));
         nt.recoElectronDzError_.push_back(track->dzError());
         nt.recoElectronTrkChi2_.push_back(track->normalizedChi2());
         nt.recoElectronTrkNumTrackerHits_.push_back(track->hitPattern().numberOfValidTrackerHits());
         nt.recoElectronTrkNumPixHits_.push_back(track->hitPattern().numberOfValidPixelHits());
         nt.recoElectronTrkNumStripHits_.push_back(track->hitPattern().numberOfValidStripHits());
      }
   }

   // Handling low-pT electrons
   std::vector<reco::GsfTrackRef> lowpt_eleTracks{};
   for (unsigned int i = 0; i < lowPtElectronHandle_->size(); i++) {
      pat::ElectronRef ele(lowPtElectronHandle_,i);
      // Cross-cleaning with regular electrons
      float mindR = 999;
      for (int k = 0; k < (int)nt.recoElectronPt_.size(); k++) {
         float dEta = ele->eta() - nt.recoElectronEta_.at(k);
         float dPhi = ele->phi() - nt.recoElectronPhi_.at(k);
         float dR = sqrt(dEta*dEta + dPhi*dPhi);
         if (dR < mindR) mindR = dR;
      }
      if (mindR < 0.01) continue;
      nt.nElectronLowPt_++;
      // Filling branches if not already in regular electron collection
      nt.recoLowPtElectronPt_.push_back(ele->pt());
      nt.recoLowPtElectronPhi_.push_back(ele->phi());
      nt.recoLowPtElectronEta_.push_back(ele->eta());
      nt.recoLowPtElectronE_.push_back(ele->energy());
      nt.recoLowPtElectronPx_.push_back(ele->px());
      nt.recoLowPtElectronPy_.push_back(ele->py());
      nt.recoLowPtElectronPz_.push_back(ele->pz());
      nt.recoLowPtElectronVxy_.push_back(ele->trackPositionAtVtx().rho());
      nt.recoLowPtElectronVz_.push_back(ele->trackPositionAtVtx().z());
      nt.recoLowPtElectronTrkIso_.push_back(ele->trackIso());
      nt.recoLowPtElectronCharge_.push_back(ele->charge());
      nt.recoLowPtElectron_passConversionVeto_.push_back(ele->passConversionVeto());
      // Filling tracks
      reco::GsfTrackRef track = ele->gsfTrack();
      if (!track.isNonnull())
         cout << "Track " << i << " from low-pT ele reco is not valid! " << endl;
      else {
         lowpt_eleTracks.emplace_back(track);
         nt.recoLowPtElectronDxy_.push_back(track->dxy(pv.position()));
         nt.recoLowPtElectronDxyError_.push_back(track->dxyError());
         nt.recoLowPtElectronDz_.push_back(track->dz(pv.position()));
         nt.recoLowPtElectronDzError_.push_back(track->dzError());
         nt.recoLowPtElectronTrkChi2_.push_back(track->normalizedChi2());
         nt.recoLowPtElectronTrkNumTrackerHits_.push_back(track->hitPattern().numberOfValidTrackerHits());
         nt.recoLowPtElectronTrkNumPixHits_.push_back(track->hitPattern().numberOfValidPixelHits());
         nt.recoLowPtElectronTrkNumStripHits_.push_back(track->hitPattern().numberOfValidStripHits());
      }
   }

   //Handling photon conversions
   for (const auto & conv : *conversionsHandle_) {
      if (conv.nTracks() != 2) continue;
      nt.nConversions_++;
      auto conv_p4 = conv.refittedPair4Momentum();
      nt.conversionPt_.push_back(conv_p4.pt());
      nt.conversionEta_.push_back(conv_p4.eta());
      nt.conversionPhi_.push_back(conv_p4.phi());
      nt.conversionE_.push_back(conv_p4.E());
      nt.conversionPx_.push_back(conv_p4.px());
      nt.conversionPy_.push_back(conv_p4.py());
      nt.conversionPz_.push_back(conv_p4.py());

      auto conv_vtx = conv.conversionVertex();
      nt.conversionVxy_.push_back(sqrt(conv_vtx.x()*conv_vtx.x() + conv_vtx.y()*conv_vtx.y()));
      nt.conversionVz_.push_back(conv_vtx.z());
      nt.conversionX_.push_back(conv_vtx.x());
      nt.conversionY_.push_back(conv_vtx.y());
      nt.conversionZ_.push_back(conv_vtx.z());

      auto t1_in = conv.tracksPin().at(0);
      auto t1_out = conv.tracksPout().at(0);
      auto t2_in = conv.tracksPin().at(1);
      auto t2_out = conv.tracksPout().at(1);
      
      nt.conversion_Trk1_innerPt_.push_back(t1_in.rho());
      nt.conversion_Trk1_innerEta_.push_back(t1_in.eta());
      nt.conversion_Trk1_innerPhi_.push_back(t1_in.phi());
      nt.conversion_Trk2_innerPt_.push_back(t2_in.rho());
      nt.conversion_Trk2_innerEta_.push_back(t2_in.eta());
      nt.conversion_Trk2_innerPhi_.push_back(t2_in.phi());
      
      nt.conversion_Trk1_outerPt_.push_back(t1_out.rho());
      nt.conversion_Trk1_outerEta_.push_back(t1_out.eta());
      nt.conversion_Trk1_outerPhi_.push_back(t1_out.phi());
      nt.conversion_Trk2_outerPt_.push_back(t2_out.rho());
      nt.conversion_Trk2_outerEta_.push_back(t2_out.eta());
      nt.conversion_Trk2_outerPhi_.push_back(t2_out.phi());
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
   auto computeVertices = [&](vector<reco::GsfTrackRef> coll_1, vector<reco::GsfTrackRef> coll_2, std::string type) {
        for (size_t i = 0; i < coll_1.size(); i++) {
            for (size_t j = 0; j < coll_2.size(); j++) {
               reco::GsfTrackRef ele_i, ele_j;
               ele_i = coll_1[i];
               ele_j = coll_2[j];
               bool skip = false;
               if ( (type == "regreg" || type == "lowlow") && i == j ) skip = true; // don't vertex ele with itself
               if (ele_i == ele_j) skip = true; // if same ele is in reg and low-pT collections

               TransientVertex tv;
               if (ele_i.isNonnull() && ele_j.isNonnull() && !skip) {
                  vector<reco::TransientTrack> transient_tracks{};
                  transient_tracks.push_back(theB->build(ele_i));
                  transient_tracks.push_back(theB->build(ele_j));
                  tv = kvf.vertex(transient_tracks);
               }
               float vxy = -9999;
               float x = -9999; float y = -9999; float z = -9999;
               float sigma_vxy = -9999;
               float vtx_chi2 = -9999;
               float vz = -9999;
               float dr = -9999;
               if (tv.isValid()) {
                  reco::Vertex vertex = reco::Vertex(tv);
                  x = vertex.x(); y = vertex.y(); z = vertex.z();
                  vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
                  sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                           vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
                  //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
                  vtx_chi2 = vertex.normalizedChi2();
                  vz = vertex.z();
                  dr = reco::deltaR(*ele_i, *ele_j);
               }

               if (type == "lowlow") {
                  pat::ElectronRef pat_ele_i(lowPtElectronHandle_,i);
                  pat::ElectronRef pat_ele_j(lowPtElectronHandle_,j);
                  auto ll = pat_ele_i->p4() + pat_ele_j->p4();
                  nt.lowlow_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.lowlow_recoVtxVxy_.push_back(vxy);
                  nt.lowlow_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.lowlow_recoVtxVz_.push_back(vz);
                  nt.lowlow_recoVtx_x_.push_back(x);
                  nt.lowlow_recoVtx_y_.push_back(y);
                  nt.lowlow_recoVtx_z_.push_back(z);
                  nt.lowlow_recoVtxDr_.push_back(dr);
                  nt.lowlow_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
                  nt.lowlow_ll_pt_.push_back(ll.pt());
                  nt.lowlow_ll_eta_.push_back(ll.eta());
                  nt.lowlow_ll_phi_.push_back(ll.phi());
                  nt.lowlow_ll_e_.push_back(ll.e());
                  nt.lowlow_ll_m_.push_back(ll.M());
                  nt.lowlow_ll_px_.push_back(ll.px());
                  nt.lowlow_ll_py_.push_back(ll.py());
                  nt.lowlow_ll_pz_.push_back(ll.pz());
               }
               else if (type == "regreg") {
                  pat::ElectronRef pat_ele_i(recoElectronHandle_,i);
                  pat::ElectronRef pat_ele_j(recoElectronHandle_,j);
                  auto ll = pat_ele_i->p4() + pat_ele_j->p4();
                  nt.regreg_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.regreg_recoVtxVxy_.push_back(vxy);
                  nt.regreg_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.regreg_recoVtxVz_.push_back(vz);
                  nt.regreg_recoVtx_x_.push_back(x);
                  nt.regreg_recoVtx_y_.push_back(y);
                  nt.regreg_recoVtx_z_.push_back(z);
                  nt.regreg_recoVtxDr_.push_back(dr);
                  nt.regreg_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
                  nt.regreg_ll_pt_.push_back(ll.pt());
                  nt.regreg_ll_eta_.push_back(ll.eta());
                  nt.regreg_ll_phi_.push_back(ll.phi());
                  nt.regreg_ll_e_.push_back(ll.e());
                  nt.regreg_ll_m_.push_back(ll.M());
                  nt.regreg_ll_px_.push_back(ll.px());
                  nt.regreg_ll_py_.push_back(ll.py());
                  nt.regreg_ll_pz_.push_back(ll.pz());
               }
               else if (type == "lowreg") {
                  pat::ElectronRef pat_ele_i(lowPtElectronHandle_,i);
                  pat::ElectronRef pat_ele_j(recoElectronHandle_,j);
                  auto ll = pat_ele_i->p4() + pat_ele_j->p4();
                  nt.lowreg_recoVtxReducedChi2_.push_back(vtx_chi2);
                  nt.lowreg_recoVtxVxy_.push_back(vxy);
                  nt.lowreg_recoVtxSigmaVxy_.push_back(sigma_vxy);
                  nt.lowreg_recoVtxVz_.push_back(vz);
                  nt.lowreg_recoVtx_x_.push_back(x);
                  nt.lowreg_recoVtx_y_.push_back(y);
                  nt.lowreg_recoVtx_z_.push_back(z);
                  nt.lowreg_recoVtxDr_.push_back(dr);
                  nt.lowreg_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
                  nt.lowreg_ll_pt_.push_back(ll.pt());
                  nt.lowreg_ll_eta_.push_back(ll.eta());
                  nt.lowreg_ll_phi_.push_back(ll.phi());
                  nt.lowreg_ll_e_.push_back(ll.e());
                  nt.lowreg_ll_m_.push_back(ll.M());
                  nt.lowreg_ll_px_.push_back(ll.px());
                  nt.lowreg_ll_py_.push_back(ll.py());
                  nt.lowreg_ll_pz_.push_back(ll.pz());
               }
            }
        }
    };

    // lowpT-lowpT
    computeVertices(lowpt_eleTracks, lowpt_eleTracks, "lowlow");
    nt.nEleVertex_lowlow_ = nt.lowlow_recoVtxVxy_.size();
    // regular-regular
    computeVertices(reg_eleTracks, reg_eleTracks, "regreg");
    nt.nEleVertex_regreg_ = nt.regreg_recoVtxVxy_.size();
    // lowpT-regular
    computeVertices(lowpt_eleTracks, reg_eleTracks, "lowreg");
    nt.nEleVertex_lowreg_ = nt.lowreg_recoVtxVxy_.size();

   ///////////////////////////////////////

   //Handling gen particles
   
   if (!isData) {
      // Gen weight
      nt.genwgt_ = genEvtInfoHandle_->weight();
      for (size_t i = 0; i < genParticleHandle_->size(); i++) {
         reco::GenParticleRef genParticle(genParticleHandle_, i);
         if (!genParticle->isHardProcess()) continue;
         nt.nGen_++;
         nt.genID_.push_back(genParticle->pdgId());
         nt.genCharge_.push_back(genParticle->charge());
         nt.genPt_.push_back(genParticle->pt());
         nt.genEta_.push_back(genParticle->eta());
         nt.genPhi_.push_back(genParticle->phi());
         nt.genEn_.push_back(genParticle->energy());
         nt.genPx_.push_back(genParticle->px());
         nt.genPy_.push_back(genParticle->py());
         nt.genPz_.push_back(genParticle->pz());
         nt.genVxy_.push_back(genParticle->vertex().rho());
         nt.genVtx_x_.push_back(genParticle->vertex().x());
         nt.genVtx_y_.push_back(genParticle->vertex().y());
         nt.genVtx_z_.push_back(genParticle->vertex().z());
         nt.genVz_.push_back(genParticle->vz());
         nt.genMass_.push_back(genParticle->mass());
      }
      
      // all gen jets
      nt.nGenJet_ = (int)genJetHandle_->size();
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
         nt.genLeadMETET_ = metRef->sumEt();
         nt.genLeadMETPx_ = metRef->px();
         nt.genLeadMETPy_ = metRef->py();
      }
    }

   outT->Fill();
   return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(iDMeSkimmer);

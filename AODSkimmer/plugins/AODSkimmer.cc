// -*- C++ -*-
//
// Package:    iDMeAnalysis/AODSkimmer
// Class:      AODSkimmer
//
/**\class AODSkimmer AODSkimmer.cc iDMeAnalysis/AODSkimmer/plugins/AODSkimmer.cc

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "iDMeAnalysis/CustomTools/interface/DisplacedDileptonAOD.hh"

#include "TTree.h"

#include "NtupleContainer.hh"

class AODSkimmer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AODSkimmer(const edm::ParameterSet&);
      ~AODSkimmer();

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
      const edm::EDGetTokenT<vector<reco::GsfElectron> > recoElectronToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> lowPtElectronToken_;
      const edm::EDGetTokenT<vector<pat::IsolatedTrack> > isoTracksToken_;
      const edm::EDGetTokenT<vector<pat::PackedCandidate> > packedPFCandToken_;
      const edm::EDGetTokenT<vector<reco::GenParticle> > genParticleToken_;
      const edm::EDGetTokenT<vector<reco::GenJet> > genJetToken_;
      const edm::EDGetTokenT<vector<reco::GenMET> > genMETToken_;
      const edm::EDGetTokenT<GenEventInfoProduct   > genEvtInfoToken_;
      const edm::EDGetTokenT<vector<reco::Vertex> > primaryVertexToken_;
      const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
      const edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
      const edm::EDGetTokenT<vector<reco::Photon> > photonsToken_;
      const edm::EDGetTokenT<vector<reco::Photon> > ootPhotonsToken_;
      const edm::EDGetTokenT<vector<reco::PFMET> > PFMETToken_;
      const edm::EDGetTokenT<vector<reco::CaloMET> > CaloMETToken_;

      // Handles
      edm::Handle<vector<reco::GsfElectron> > recoElectronHandle_;
      edm::Handle<pat::ElectronCollection> lowPtElectronHandle_;
      edm::Handle<vector<pat::IsolatedTrack> > isoTracksHandle_;
      edm::Handle<vector<pat::PackedCandidate> > packedPFCandHandle_;
      edm::Handle<vector<reco::GenParticle> > genParticleHandle_;
      edm::Handle<vector<reco::GenJet> > genJetHandle_;
      edm::Handle<vector<reco::GenMET> > genMETHandle_;
      edm::Handle<GenEventInfoProduct> genEvtInfoHandle_;
      edm::Handle<vector<reco::Vertex> > primaryVertexHandle_;
      edm::Handle<reco::BeamSpot> beamspotHandle_;
      edm::Handle<vector<reco::Conversion> > conversionsHandle_;
      edm::Handle<vector<reco::Photon> > photonsHandle_;
      edm::Handle<vector<reco::Photon> > ootPhotonsHandle_;
      edm::Handle<vector<reco::PFMET> > PFMETHandle_;
      edm::Handle<vector<reco::CaloMET> > CaloMETHandle_;
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
AODSkimmer::AODSkimmer(const edm::ParameterSet& ps)
 :
   isData(ps.getParameter<bool>("isData")),
   recoElectronToken_(consumes<vector<reco::GsfElectron> >(ps.getParameter<edm::InputTag>("recoElectron"))),
   lowPtElectronToken_(consumes<pat::ElectronCollection>(ps.getParameter<edm::InputTag>("lowPtElectron"))),
   isoTracksToken_(consumes<vector<pat::IsolatedTrack> >(ps.getParameter<edm::InputTag>("isoTracks"))),
   packedPFCandToken_(consumes<vector<pat::PackedCandidate> >(ps.getParameter<edm::InputTag>("packedPFCands"))),
   genParticleToken_(consumes<vector<reco::GenParticle> >(ps.getParameter<edm::InputTag>("genParticle"))),
   genJetToken_(consumes<vector<reco::GenJet> >(ps.getParameter<edm::InputTag>("genJet"))),
   genMETToken_(consumes<vector<reco::GenMET> >(ps.getParameter<edm::InputTag>("genMET"))),
   genEvtInfoToken_(consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("genEvt"))),
   primaryVertexToken_(consumes<vector<reco::Vertex> >(ps.getParameter<edm::InputTag>("primaryVertex"))),
   beamspotToken_(consumes<reco::BeamSpot>(ps.getParameter<edm::InputTag>("beamSpot"))),
   conversionsToken_(consumes<vector<reco::Conversion> >(ps.getParameter<edm::InputTag>("conversions"))),
   photonsToken_(consumes<vector<reco::Photon> >(ps.getParameter<edm::InputTag>("photons"))),
   ootPhotonsToken_(consumes<vector<reco::Photon> >(ps.getParameter<edm::InputTag>("ootPhotons"))),
   PFMETToken_(consumes<vector<reco::PFMET> >(ps.getParameter<edm::InputTag>("PFMET"))),
   CaloMETToken_(consumes<vector<reco::CaloMET> >(ps.getParameter<edm::InputTag>("CaloMET")))
{
   usesResource("TFileService");
   m_random_generator = std::mt19937(37428479);

}


AODSkimmer::~AODSkimmer() = default;


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void
AODSkimmer::beginJob()
{
   outT = fs->make<TTree>("outT", "outT");
   nt.isData_ = isData;
   nt.SetTree(outT);
   nt.CreateTreeBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void
AODSkimmer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AODSkimmer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   edm::ParameterSetDescription desc;

   desc.add<bool>("isData", 0);

   desc.add<edm::InputTag>("recoElectron",edm::InputTag("gedGsfElectrons"));
   desc.add<edm::InputTag>("lowPtElectron",edm::InputTag("slimmedLowPtElectrons"));
   desc.add<edm::InputTag>("isoTracks",edm::InputTag("isolatedTracks"));
   desc.add<edm::InputTag>("packedPFCands",edm::InputTag("packedPFCandidates"));
   desc.add<edm::InputTag>("genParticle",edm::InputTag("genParticles"));
   desc.add<edm::InputTag>("genJet",edm::InputTag("ak4GenJets"));
   desc.add<edm::InputTag>("genMET",edm::InputTag("genMetTrue"));
   desc.add<edm::InputTag>("genEvt", edm::InputTag("generator"));
   desc.add<edm::InputTag>("primaryVertex",edm::InputTag("offlinePrimaryVertices"));
   desc.add<edm::InputTag>("beamSpot",edm::InputTag("offlineBeamSpot"));
   desc.add<edm::InputTag>("conversions",edm::InputTag("allConversions"));
   desc.add<edm::InputTag>("photons",edm::InputTag("photons"));
   desc.add<edm::InputTag>("ootPhotons",edm::InputTag("ootPhotons"));
   desc.add<edm::InputTag>("PFMET",edm::InputTag("pfMet"));
   desc.add<edm::InputTag>("CaloMET",edm::InputTag("caloMet"));
   
   descriptions.add("AODSkimmer", desc);
}

// ------------ method called for each event  ------------
void
AODSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using std::cout, std::endl; 

   // Retrieving event data and assigning to handles
   iEvent.getByToken(recoElectronToken_,recoElectronHandle_);
   iEvent.getByToken(lowPtElectronToken_,lowPtElectronHandle_);
   iEvent.getByToken(isoTracksToken_,isoTracksHandle_);
   iEvent.getByToken(packedPFCandToken_,packedPFCandHandle_);
   iEvent.getByToken(genParticleToken_,genParticleHandle_);
   iEvent.getByToken(genJetToken_,genJetHandle_);
   iEvent.getByToken(genMETToken_,genMETHandle_);
   iEvent.getByToken(genEvtInfoToken_,genEvtInfoHandle_);
   iEvent.getByToken(primaryVertexToken_,primaryVertexHandle_);
   iEvent.getByToken(beamspotToken_,beamspotHandle_);
   iEvent.getByToken(conversionsToken_,conversionsHandle_);
   iEvent.getByToken(photonsToken_,photonsHandle_);
   iEvent.getByToken(ootPhotonsToken_,ootPhotonsHandle_);
   iEvent.getByToken(PFMETToken_,PFMETHandle_);
   iEvent.getByToken(CaloMETToken_,CaloMETHandle_);

   // Clear tree branches before filling
   nt.ClearTreeBranches();

   //////
   // Computing derived quantities and filling the trees
   //////
   nt.eventNum_ = iEvent.id().event();
   nt.lumiSec_ = iEvent.luminosityBlock();
   nt.runNum_ = iEvent.id().run();
   reco::Vertex pv = (*primaryVertexHandle_).at(0);
   reco::BeamSpot beamspot = *beamspotHandle_;
   // Set up objects for vertex reco
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   KalmanVertexFitter kvf(true);

   // Handling Regular MET
   
   if (PFMETHandle_->size() > 0) {
      auto pfmet = (*PFMETHandle_).at(0);
      nt.PFMET_ET_ = pfmet.sumEt();
      nt.PFMET_Px_ = pfmet.px();
      nt.PFMET_Py_ = pfmet.py();
      nt.PFMET_Pt_ = pfmet.pt();
      nt.PFMET_Phi_ = pfmet.phi();
   }

   if (CaloMETHandle_->size() > 0) {
      auto calomet = (*CaloMETHandle_).at(0);
      nt.CaloMET_ET_ = calomet.sumEt();
      nt.CaloMET_Px_ = calomet.px();
      nt.CaloMET_Py_ = calomet.py();
      nt.CaloMET_Pt_ = calomet.pt();
      nt.CaloMET_Phi_ = calomet.phi();
   }
   
   // Handling default electrons
   nt.nElectronDefault_ = recoElectronHandle_->size();
   vector<reco::GsfTrackRef> reg_eleTracks{};
   vector<math::XYZTLorentzVector> reg_ele_p4s;
   for (const auto & ele : *recoElectronHandle_) {
      nt.recoElectronPt_.push_back(ele.pt());
      nt.recoElectronEta_.push_back(ele.eta());
      nt.recoElectronPhi_.push_back(ele.phi());
      nt.recoElectronE_.push_back(ele.energy());
      nt.recoElectronPx_.push_back(ele.px());
      nt.recoElectronPy_.push_back(ele.py());
      nt.recoElectronPz_.push_back(ele.pz());
      nt.recoElectronVxy_.push_back(ele.trackPositionAtVtx().rho());
      nt.recoElectronVz_.push_back(ele.trackPositionAtVtx().z());
      nt.recoElectronTrkIso_.push_back(ele.dr04TkSumPt());
      nt.recoElectronCharge_.push_back(ele.charge());
      // Filling tracks
      reco::GsfTrackRef track = ele.gsfTrack();
      reg_eleTracks.push_back(track);
      reg_ele_p4s.push_back(ele.p4());
      nt.recoElectronDxy_.push_back(track->dxy(pv.position()));
      nt.recoElectronDxyError_.push_back(track->dxyError());
      nt.recoElectronDz_.push_back(track->dz(pv.position()));
      nt.recoElectronDzError_.push_back(track->dzError());
      nt.recoElectronTrkChi2_.push_back(track->normalizedChi2());
      nt.recoElectronTrkNumTrackerHits_.push_back(track->hitPattern().numberOfValidTrackerHits());
      nt.recoElectronTrkNumPixHits_.push_back(track->hitPattern().numberOfValidPixelHits());
      nt.recoElectronTrkNumStripHits_.push_back(track->hitPattern().numberOfValidStripHits());
   }

   // Handling low-pT electrons
   std::vector<reco::GsfTrackRef> lowpt_eleTracks{};
   vector<math::XYZTLorentzVector> lowpt_ele_p4s;
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
      lowpt_eleTracks.push_back(track);
      lowpt_ele_p4s.push_back(ele->p4());
      nt.recoLowPtElectronDxy_.push_back(track->dxy(pv.position()));
      nt.recoLowPtElectronDxyError_.push_back(track->dxyError());
      nt.recoLowPtElectronDz_.push_back(track->dz(pv.position()));
      nt.recoLowPtElectronDzError_.push_back(track->dzError());
      nt.recoLowPtElectronTrkChi2_.push_back(track->normalizedChi2());
      nt.recoLowPtElectronTrkNumTrackerHits_.push_back(track->hitPattern().numberOfValidTrackerHits());
      nt.recoLowPtElectronTrkNumPixHits_.push_back(track->hitPattern().numberOfValidPixelHits());
      nt.recoLowPtElectronTrkNumStripHits_.push_back(track->hitPattern().numberOfValidStripHits());
   }

   // Handling photons
   for (const auto & ph : *photonsHandle_) {
      nt.nPhotons_++;
      nt.PhotonEt_.push_back(ph.et());
      nt.PhotonEt_.push_back(ph.eta());
      nt.PhotonEt_.push_back(ph.phi());
   }

   // Handling OOT photons
   for (const auto & ph : *ootPhotonsHandle_) {
      nt.nOOTPhotons_++;
      nt.ootPhotonEt_.push_back(ph.et());
      nt.ootPhotonEt_.push_back(ph.eta());
      nt.ootPhotonEt_.push_back(ph.phi());
   }

   //Handling photon conversions
   for (const auto & conv : *conversionsHandle_) {
      if (conv.nTracks() == 0) continue;
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
      nt.conversion_Trk1_innerPt_.push_back(t1_in.rho());
      nt.conversion_Trk1_innerEta_.push_back(t1_in.eta());
      nt.conversion_Trk1_innerPhi_.push_back(t1_in.phi());
      nt.conversion_Trk1_outerPt_.push_back(t1_out.rho());
      nt.conversion_Trk1_outerEta_.push_back(t1_out.eta());
      nt.conversion_Trk1_outerPhi_.push_back(t1_out.phi());
      
      if (conv.nTracks() > 1) {
         auto t2_in = conv.tracksPin().at(1);
         auto t2_out = conv.tracksPout().at(1);
         nt.conversion_Trk2_innerPt_.push_back(t2_in.rho());
         nt.conversion_Trk2_innerEta_.push_back(t2_in.eta());
         nt.conversion_Trk2_innerPhi_.push_back(t2_in.phi());
         nt.conversion_Trk2_outerPt_.push_back(t2_out.rho());
         nt.conversion_Trk2_outerEta_.push_back(t2_out.eta());
         nt.conversion_Trk2_outerPhi_.push_back(t2_out.phi());
      }
      else {
         nt.conversion_Trk2_innerPt_.push_back(-999);
         nt.conversion_Trk2_innerEta_.push_back(-999);
         nt.conversion_Trk2_innerPhi_.push_back(-999);
         nt.conversion_Trk2_outerPt_.push_back(-999);
         nt.conversion_Trk2_outerEta_.push_back(-999);
         nt.conversion_Trk2_outerPhi_.push_back(-999);
      }
   }

   // Algorithm to reconstruct displaced dilepton vertices from isoTracks
   // Want to use regular and OOT photons for this
   vector<reco::Photon> allPhotons = *photonsHandle_;
   allPhotons.insert(allPhotons.end(),(*ootPhotonsHandle_).begin(),(*ootPhotonsHandle_).end());
   DisplacedDileptonsAOD EleFinder(pv,beamspot,isoTracksHandle_,allPhotons,packedPFCandHandle_,theB);
   EleFinder.findDileptons();
   
   // Extracting tracks from new electron candidates
   vector<reco::Track> cand_eleTracks;
   vector<math::XYZTLorentzVector> cand_ele_p4s;
   for (int icand = 0; icand < EleFinder.nElectronCandidate; icand++) {
      int tk_idx = EleFinder.ElectronCandidate_isotrackIdx[icand];
      auto isoTk = (*isoTracksHandle_).at(tk_idx);
      auto tk = *(isoTk.packedCandRef()->bestTrack());
      cand_eleTracks.push_back(tk);
      cand_ele_p4s.push_back(isoTk.p4());
   }
   
   // Filling additional electron candidates
   nt.nEleCand_ = EleFinder.nElectronCandidate;
   for (int iec = 0; iec < nt.nEleCand_; iec++) {
      nt.EleCand_pt_.push_back(EleFinder.ElectronCandidate_pt[iec]);
      nt.EleCand_et_.push_back(EleFinder.ElectronCandidate_et[iec]);
      nt.EleCand_eta_.push_back(EleFinder.ElectronCandidate_eta[iec]);
      nt.EleCand_phi_.push_back(EleFinder.ElectronCandidate_phi[iec]);
      nt.EleCand_dxy_.push_back(EleFinder.ElectronCandidate_dxy[iec]);
      nt.EleCand_dxyError_.push_back(EleFinder.ElectronCandidate_dxyError[iec]);
      nt.EleCand_dxy_PV_.push_back(EleFinder.ElectronCandidate_dxy_PV[iec]);
      nt.EleCand_dxyError_PV_.push_back(EleFinder.ElectronCandidate_dxyError_PV[iec]);
      nt.EleCand_relPFiso_.push_back(EleFinder.ElectronCandidate_relPFiso[iec]);
      nt.EleCand_relTrkiso_.push_back(EleFinder.ElectronCandidate_relTrkiso[iec]);
      nt.EleCand_ptDiff_.push_back(EleFinder.ElectronCandidate_ptDiff[iec]);
      nt.EleCand_trkIso_.push_back(EleFinder.ElectronCandidate_trkIso[iec]);
      nt.EleCand_trkChi2_.push_back(EleFinder.ElectronCandidate_trkChi2[iec]);
      nt.EleCand_numTrackerHits_.push_back(EleFinder.ElectronCandidate_numTrackerHits[iec]);
      nt.EleCand_numPixHits_.push_back(EleFinder.ElectronCandidate_numPixHits[iec]);
      nt.EleCand_numStripHits_.push_back(EleFinder.ElectronCandidate_numStripHits[iec]);
   }
   
   // Filling displaced dileptons that pass baseline selection
   nt.ndispEE_ = EleFinder.nEEBase;
   nt.dispEE_maxIxy_.push_back(EleFinder.EEBase_maxIxy);
   for (int iee = 0; iee < nt.ndispEE_; iee++) {
      nt.dispEE_Lxy_.push_back(EleFinder.EEBase_Lxy[iee]);
      nt.dispEE_Ixy_.push_back(EleFinder.EEBase_Ixy[iee]);
      nt.dispEE_trackDxy_.push_back(EleFinder.EEBase_trackDxy[iee]);
      nt.dispEE_trackIxy_.push_back(EleFinder.EEBase_trackIxy[iee]);
      nt.dispEE_vx_.push_back(EleFinder.EEBase_vx[iee]);
      nt.dispEE_vy_.push_back(EleFinder.EEBase_vy[iee]);
      nt.dispEE_mass_.push_back(EleFinder.EEBase_mass[iee]);
      nt.dispEE_normalizedChi2_.push_back(EleFinder.EEBase_normalizedChi2[iee]);
      nt.dispEE_leadingPt_.push_back(EleFinder.EEBase_leadingPt[iee]);
      nt.dispEE_subleadingPt_.push_back(EleFinder.EEBase_subleadingPt[iee]);
      nt.dispEE_leadingEt_.push_back(EleFinder.EEBase_leadingEt[iee]);
      nt.dispEE_subleadingEt_.push_back(EleFinder.EEBase_subleadingEt[iee]);
      nt.dispEE_cosAlpha_.push_back(EleFinder.EEBase_cosAlpha[iee]);
      nt.dispEE_dPhi_.push_back(EleFinder.EEBase_dPhi[iee]);
      nt.dispEE_relisoA_.push_back(EleFinder.EEBase_relisoA[iee]);
      nt.dispEE_relisoB_.push_back(EleFinder.EEBase_relisoB[iee]);
      nt.dispEE_fromPVA_.push_back(EleFinder.EEBase_fromPVA[iee]);
      nt.dispEE_fromPVB_.push_back(EleFinder.EEBase_fromPVB[iee]);
      nt.dispEE_PVAssociation_.push_back(EleFinder.EEBase_PVAssociation[iee]);
   }

   // Filling displaced dilepton candidates
   nt.nEECand_ = EleFinder.nEE;
   for (int ieec = 0; ieec < nt.nEECand_; ieec++) {
      nt.EECand_Lxy_PV_.push_back(EleFinder.EE_Lxy_PV[ieec]);
      nt.EECand_Ixy_PV_.push_back(EleFinder.EE_Ixy_PV[ieec]);
      nt.EECand_Lxy_0_.push_back(EleFinder.EE_Lxy_0[ieec]);
      nt.EECand_Ixy_0_.push_back(EleFinder.EE_Ixy_0[ieec]);
      nt.EECand_Lxy_BS_.push_back(EleFinder.EE_Lxy_BS[ieec]);
      nt.EECand_Ixy_BS_.push_back(EleFinder.EE_Ixy_BS[ieec]);
      nt.EECand_trackDxy_.push_back(EleFinder.EE_trackDxy[ieec]);
      nt.EECand_trackIxy_.push_back(EleFinder.EE_trackIxy[ieec]);
      nt.EECand_trackDxy_PV_.push_back(EleFinder.EE_trackDxy_PV[ieec]);
      nt.EECand_trackIxy_PV_.push_back(EleFinder.EE_trackIxy_PV[ieec]);
      nt.EECand_trackDxy_0_.push_back(EleFinder.EE_trackDxy_0[ieec]);
      nt.EECand_trackIxy_0_.push_back(EleFinder.EE_trackIxy_0[ieec]);
      nt.EECand_trackDxy_BS_.push_back(EleFinder.EE_trackDxy_BS[ieec]);
      nt.EECand_trackIxy_BS_.push_back(EleFinder.EE_trackIxy_BS[ieec]);
      nt.EECand_vx_.push_back(EleFinder.EE_vx[ieec]);
      nt.EECand_vy_.push_back(EleFinder.EE_vy[ieec]);
      nt.EECand_mass_.push_back(EleFinder.EE_mass[ieec]);
      nt.EECand_normalizedChi2_.push_back(EleFinder.EE_normalizedChi2[ieec]);
      nt.EECand_leadingPt_.push_back(EleFinder.EE_leadingPt[ieec]);
      nt.EECand_subleadingPt_.push_back(EleFinder.EE_subleadingPt[ieec]);
      nt.EECand_leadingEt_.push_back(EleFinder.EE_leadingEt[ieec]);
      nt.EECand_subleadingEt_.push_back(EleFinder.EE_subleadingEt[ieec]);
      nt.EECand_cosAlpha_.push_back(EleFinder.EE_cosAlpha[ieec]);
      nt.EECand_dR_.push_back(EleFinder.EE_dR[ieec]);
      nt.EECand_dPhi_.push_back(EleFinder.EE_dPhi[ieec]);
      nt.EECand_lldPhi_.push_back(EleFinder.EE_lldPhi[ieec]);
      nt.EECand_relisoA_.push_back(EleFinder.EE_relisoA[ieec]);
      nt.EECand_relisoB_.push_back(EleFinder.EE_relisoB[ieec]);
   }

   // Define vertex reco function 
   auto computeVertices = [&](vector<int> coll_1, vector<int> coll_2, std::string type) {
      for (size_t i = 0; i < coll_1.size(); i++) {
         for (size_t j = 0; j < coll_2.size(); j++) {
            int ind1 = coll_1[i];
            int ind2 = coll_2[j];
            reco::GsfTrackRef ele_i, ele_j;
            if (type == "regreg") {
               ele_i = reg_eleTracks[ind1];
               ele_j = reg_eleTracks[ind2];
            }
            else if (type == "lowlow") {
               ele_i = lowpt_eleTracks[ind1];
               ele_j = lowpt_eleTracks[ind2];
            }
            else if (type == "lowreg") {
               ele_i = lowpt_eleTracks[ind1];
               ele_j = reg_eleTracks[ind2];
            }
            if ( (type == "regreg" || type == "lowlow") && (j <= i) ) continue; // don't vertex ele with itself or ones prior (if vertexing with same type)
            if (ele_i == ele_j) continue; // skip if same ele is in reg and low-pT collections
            if (!ele_i.isNonnull() || !ele_j.isNonnull()) continue; // skip if there's a bad track

            TransientVertex tv;
            vector<reco::TransientTrack> transient_tracks{};
            transient_tracks.push_back(theB->build(ele_i));
            transient_tracks.push_back(theB->build(ele_j));
            tv = kvf.vertex(transient_tracks);

            if (!tv.isValid()) continue; // skip if the vertex is bad

            reco::Vertex vertex = reco::Vertex(tv);
            float vx = vertex.x(); 
            float vy = vertex.y(); 
            float vz = vertex.z();
            float vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
            float sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                     vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
            float vtx_chi2 = vertex.normalizedChi2();
            float dr = reco::deltaR(*ele_i, *ele_j);

            if (type == "lowlow") {
               math::XYZTLorentzVector ll = lowpt_ele_p4s[ind1] + lowpt_ele_p4s[ind2];
               nt.LLvtx_idx1_.push_back(ind1); 
               nt.LLvtx_idx2_.push_back(ind2);
               nt.LLvtx_recoVtxReducedChi2_.push_back(vtx_chi2);
               nt.LLvtx_recoVtxVxy_.push_back(vxy);
               nt.LLvtx_recoVtxSigmaVxy_.push_back(sigma_vxy);
               nt.LLvtx_recoVtxVx_.push_back(vx);
               nt.LLvtx_recoVtxVy_.push_back(vy);
               nt.LLvtx_recoVtxVz_.push_back(vz);
               nt.LLvtx_recoVtxDr_.push_back(dr);
               nt.LLvtx_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
               nt.LLvtx_ll_pt_.push_back(ll.pt());
               nt.LLvtx_ll_eta_.push_back(ll.eta());
               nt.LLvtx_ll_phi_.push_back(ll.phi());
               nt.LLvtx_ll_e_.push_back(ll.e());
               nt.LLvtx_ll_m_.push_back(ll.M());
               nt.LLvtx_ll_px_.push_back(ll.px());
               nt.LLvtx_ll_py_.push_back(ll.py());
               nt.LLvtx_ll_pz_.push_back(ll.pz());
            }
            else if (type == "regreg") {
               math::XYZTLorentzVector ll = reg_ele_p4s[ind1] + reg_ele_p4s[ind2];
               nt.RRvtx_idx1_.push_back(ind1); 
               nt.RRvtx_idx2_.push_back(ind2);
               nt.RRvtx_recoVtxReducedChi2_.push_back(vtx_chi2);
               nt.RRvtx_recoVtxVxy_.push_back(vxy);
               nt.RRvtx_recoVtxSigmaVxy_.push_back(sigma_vxy);
               nt.RRvtx_recoVtxVx_.push_back(vx);
               nt.RRvtx_recoVtxVy_.push_back(vy);
               nt.RRvtx_recoVtxVz_.push_back(vz);
               nt.RRvtx_recoVtxDr_.push_back(dr);
               nt.RRvtx_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
               nt.RRvtx_ll_pt_.push_back(ll.pt());
               nt.RRvtx_ll_eta_.push_back(ll.eta());
               nt.RRvtx_ll_phi_.push_back(ll.phi());
               nt.RRvtx_ll_e_.push_back(ll.e());
               nt.RRvtx_ll_m_.push_back(ll.M());
               nt.RRvtx_ll_px_.push_back(ll.px());
               nt.RRvtx_ll_py_.push_back(ll.py());
               nt.RRvtx_ll_pz_.push_back(ll.pz());
            }
            else if (type == "lowreg") {
               math::XYZTLorentzVector ll = lowpt_ele_p4s[ind1] + reg_ele_p4s[ind2];
               nt.LRvtx_idx1_.push_back(ind1); 
               nt.LRvtx_idx2_.push_back(ind2);
               nt.LRvtx_recoVtxReducedChi2_.push_back(vtx_chi2);
               nt.LRvtx_recoVtxVxy_.push_back(vxy);
               nt.LRvtx_recoVtxSigmaVxy_.push_back(sigma_vxy);
               nt.LRvtx_recoVtxVx_.push_back(vx);
               nt.LRvtx_recoVtxVy_.push_back(vy);
               nt.LRvtx_recoVtxVz_.push_back(vz);
               nt.LRvtx_recoVtxDr_.push_back(dr);
               nt.LRvtx_recoVtxSign_.push_back(ele_i->charge()*ele_j->charge());
               nt.LRvtx_ll_pt_.push_back(ll.pt());
               nt.LRvtx_ll_eta_.push_back(ll.eta());
               nt.LRvtx_ll_phi_.push_back(ll.phi());
               nt.LRvtx_ll_e_.push_back(ll.e());
               nt.LRvtx_ll_m_.push_back(ll.M());
               nt.LRvtx_ll_px_.push_back(ll.px());
               nt.LRvtx_ll_py_.push_back(ll.py());
               nt.LRvtx_ll_pz_.push_back(ll.pz());
            }
         }
      }
   };

   auto computeVerticesCand = [&](vector<int> coll_1, vector<int> coll_2, std::string type) {
      for (size_t i = 0; i < coll_1.size(); i++) {
         for (size_t j = 0; j < coll_2.size(); j++) {
            int ind1 = coll_1[i];
            int ind2 = coll_2[j];
            reco::GsfTrackRef ele_i;
            reco::Track ele_j;
            if (type == "regcand") {
               ele_i = reg_eleTracks[ind1];
               ele_j = cand_eleTracks[ind2];
            }
            else if (type == "lowcand") {
               ele_i = lowpt_eleTracks[ind1];
               ele_j = cand_eleTracks[ind2];
            }
            if (!ele_i.isNonnull()) continue; // skip if electron is bad
            if (reco::deltaR(*ele_i,ele_j) < 0.0001) continue; // skip if it's somehow the same track

            TransientVertex tv;
            vector<reco::TransientTrack> transient_tracks{};
            transient_tracks.push_back(theB->build(ele_i));
            transient_tracks.push_back(theB->build(ele_j));
            tv = kvf.vertex(transient_tracks);

            if (!tv.isValid()) continue; // skip if the vertex is bad

            reco::Vertex vertex = reco::Vertex(tv);
            float vx = vertex.x(); 
            float vy = vertex.y(); 
            float vz = vertex.z();
            float vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
            float sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
                     vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
            float vtx_chi2 = vertex.normalizedChi2();
            float dr = reco::deltaR(*ele_i, ele_j);

            if (type == "regcand") {
               math::XYZTLorentzVector ll = reg_ele_p4s[ind1] + cand_ele_p4s[ind2];
               nt.RCvtx_idx1_.push_back(ind1); 
               nt.RCvtx_idx2_.push_back(ind2);
               nt.RCvtx_recoVtxReducedChi2_.push_back(vtx_chi2);
               nt.RCvtx_recoVtxVxy_.push_back(vxy);
               nt.RCvtx_recoVtxSigmaVxy_.push_back(sigma_vxy);
               nt.RCvtx_recoVtxVx_.push_back(vx);
               nt.RCvtx_recoVtxVy_.push_back(vy);
               nt.RCvtx_recoVtxVz_.push_back(vz);
               nt.RCvtx_recoVtxDr_.push_back(dr);
               nt.RCvtx_recoVtxSign_.push_back(ele_i->charge()*ele_j.charge());
               nt.RCvtx_ll_pt_.push_back(ll.pt());
               nt.RCvtx_ll_eta_.push_back(ll.eta());
               nt.RCvtx_ll_phi_.push_back(ll.phi());
               nt.RCvtx_ll_e_.push_back(ll.e());
               nt.RCvtx_ll_m_.push_back(ll.M());
               nt.RCvtx_ll_px_.push_back(ll.px());
               nt.RCvtx_ll_py_.push_back(ll.py());
               nt.RCvtx_ll_pz_.push_back(ll.pz());
            }
            else if (type == "lowcand") {
               math::XYZTLorentzVector ll = lowpt_ele_p4s[ind1] + cand_ele_p4s[ind2];
               nt.LCvtx_idx1_.push_back(ind1); 
               nt.LCvtx_idx2_.push_back(ind2);
               nt.LCvtx_recoVtxReducedChi2_.push_back(vtx_chi2);
               nt.LCvtx_recoVtxVxy_.push_back(vxy);
               nt.LCvtx_recoVtxSigmaVxy_.push_back(sigma_vxy);
               nt.LCvtx_recoVtxVx_.push_back(vx);
               nt.LCvtx_recoVtxVy_.push_back(vy);
               nt.LCvtx_recoVtxVz_.push_back(vz);
               nt.LCvtx_recoVtxDr_.push_back(dr);
               nt.LCvtx_recoVtxSign_.push_back(ele_i->charge()*ele_j.charge());
               nt.LCvtx_ll_pt_.push_back(ll.pt());
               nt.LCvtx_ll_eta_.push_back(ll.eta());
               nt.LCvtx_ll_phi_.push_back(ll.phi());
               nt.LCvtx_ll_e_.push_back(ll.e());
               nt.LCvtx_ll_m_.push_back(ll.M());
               nt.LCvtx_ll_px_.push_back(ll.px());
               nt.LCvtx_ll_py_.push_back(ll.py());
               nt.LCvtx_ll_pz_.push_back(ll.pz());
            }
         }
      }
   };

   // Reconstructing electron vertices
   // Sort electron tracks by pT
   vector<int> ind_ele(reg_eleTracks.size());
   std::iota(ind_ele.begin(),ind_ele.end(),0); // create vector of electron indices
   std::sort(ind_ele.begin(), ind_ele.end(), [&](const int & l, const int & r) { // sort indices by ele pT
            return reg_ele_p4s[l].pt() > reg_ele_p4s[r].pt();
            });
   // Sort low-pT electron tracks by pT
   vector<int> ind_lptele(lowpt_eleTracks.size());
   std::iota(ind_lptele.begin(),ind_lptele.end(),0); // create vector of low-pT electron indices
   std::sort(ind_lptele.begin(), ind_lptele.end(), [&](const int & l, const int & r) { // sort indices by ele pT
            return lowpt_ele_p4s[l].pt() > lowpt_ele_p4s[r].pt();
            });
   // Sort isotrack electron candidates by pT
   vector<int> ind_eleCand(cand_eleTracks.size());
   std::iota(ind_eleCand.begin(),ind_eleCand.end(),0);
   std::sort(ind_eleCand.begin(),ind_eleCand.end(), [&](const int & l, const int &r) {
      return cand_ele_p4s[l].pt() > cand_ele_p4s[r].pt();
   });


   // lowpT-lowpT
   computeVertices(ind_lptele, ind_lptele, "lowlow");
   nt.nEleVertex_LL_ = nt.LLvtx_recoVtxVxy_.size();
   // regular-regular
   computeVertices(ind_ele, ind_ele, "regreg");
   nt.nEleVertex_RR_ = nt.RRvtx_recoVtxVxy_.size();
   // lowpT-regular
   computeVertices(ind_lptele, ind_ele, "lowreg");
   nt.nEleVertex_LR_ = nt.LRvtx_recoVtxVxy_.size();
   // regular-cand
   computeVerticesCand(ind_ele,ind_eleCand,"regcand");
   nt.nEleVertex_RC_ = nt.RCvtx_recoVtxVxy_.size();
   // lowpT-cand
   computeVerticesCand(ind_lptele,ind_eleCand,"lowcand");
   nt.nEleVertex_LC_ = nt.LCvtx_recoVtxVxy_.size();

   ///////////////////////////////////////

   //Handling gen particles
   if (!isData) {
      // Gen weight
      nt.genwgt_ = genEvtInfoHandle_->weight();
      for (const auto & genParticle : *genParticleHandle_) {
         if (!genParticle.isHardProcess()) continue;
         nt.nGen_++;
         nt.genID_.push_back(genParticle.pdgId());
         nt.genCharge_.push_back(genParticle.charge());
         nt.genPt_.push_back(genParticle.pt());
         nt.genEta_.push_back(genParticle.eta());
         nt.genPhi_.push_back(genParticle.phi());
         nt.genEn_.push_back(genParticle.energy());
         nt.genPx_.push_back(genParticle.px());
         nt.genPy_.push_back(genParticle.py());
         nt.genPz_.push_back(genParticle.pz());
         nt.genVxy_.push_back(sqrt(genParticle.vx()*genParticle.vx() + genParticle.vy()*genParticle.vy()));
         nt.genVx_.push_back(genParticle.vx());
         nt.genVy_.push_back(genParticle.vy());
         nt.genVz_.push_back(genParticle.vz());
         nt.genMass_.push_back(genParticle.mass());
      }
      
      // all gen jets
      nt.nGenJet_ = (int)genJetHandle_->size();
      for (const auto & jet : *genJetHandle_) {
         nt.genJetPt_.push_back(jet.pt());
         nt.genJetEta_.push_back(jet.eta());
         nt.genJetPhi_.push_back(jet.phi());
      }

      // Lead gen MET
      if (genMETHandle_->size() > 0) {
         auto met = (*genMETHandle_).at(0);
         nt.genLeadMETPt_ = met.pt();
         nt.genLeadMETPhi_ = met.phi();
         nt.genLeadMETET_ = met.sumEt();
         nt.genLeadMETPx_ = met.px();
         nt.genLeadMETPy_ = met.py();
      }
    }

   outT->Fill();
   return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODSkimmer);

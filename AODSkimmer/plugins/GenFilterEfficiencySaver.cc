// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "CommonTools/UtilAlgos/interface/TFileService.h"
 #include <TTree.h>

#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"

#include <iostream>

class GenFilterEfficiencySaver : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchLuminosityBlocks> {
  public:
    explicit GenFilterEfficiencySaver(const edm::ParameterSet&);
    ~GenFilterEfficiencySaver();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
  private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    void endJob() override;

    edm::EDGetTokenT<GenFilterInfo> genFilterInfoToken_;
    GenFilterInfo totalGenFilterInfo_;
    TTree *outT;
    float total_evts;
    float passed_evts;
    float filter_eff;
    float filter_eff_uncert;

    // ----------member data ---------------------------
  
};

GenFilterEfficiencySaver::GenFilterEfficiencySaver(const edm::ParameterSet& pset):
  genFilterInfoToken_(consumes<GenFilterInfo,edm::InLumi>(pset.getParameter<edm::InputTag>("genFilterInfoTag"))),
  totalGenFilterInfo_(0,0,0,0,0.,0.,0.,0.)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  outT = fs->make<TTree>("outT", "outT");
  outT->Branch("total_evts",&total_evts);
  outT->Branch("passed_evts",&passed_evts);
  outT->Branch("filter_eff",&filter_eff);
  outT->Branch("fitler_eff_uncert",&filter_eff_uncert);
}

GenFilterEfficiencySaver::~GenFilterEfficiencySaver()
{
}

void GenFilterEfficiencySaver::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genFilterInfoTag",edm::InputTag("genFilterEfficiencyProducer"));
  descriptions.add("GenFilterEfficiencySaver", desc);
}

void
GenFilterEfficiencySaver::analyze(const edm::Event&, const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------

void
GenFilterEfficiencySaver::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

  edm::Handle<GenFilterInfo> genFilter;
  iLumi.getByToken(genFilterInfoToken_, genFilter);

  std::cout << "Lumi section " << iLumi.id() << std::endl;

  std::cout << "N total = " << genFilter->sumWeights() << " N passed = " << genFilter->sumPassWeights() << " N failed = " << genFilter->sumFailWeights() << std::endl;
  std::cout << "Generator filter efficiency = " << genFilter->filterEfficiency(-1) << " +- " << genFilter->filterEfficiencyError(-1) << std::endl;
  totalGenFilterInfo_.mergeProduct(*genFilter);
  
}

void
GenFilterEfficiencySaver::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void
GenFilterEfficiencySaver::endJob() {

  std::cout << "Total events = " << totalGenFilterInfo_.sumWeights()
	    << " Passed events = " << totalGenFilterInfo_.sumPassWeights() << std::endl;
  std::cout << "Filter efficiency = " << totalGenFilterInfo_.filterEfficiency(-1)  
	    << " +- " << totalGenFilterInfo_.filterEfficiencyError(-1)  << std::endl;

  total_evts = totalGenFilterInfo_.sumWeights();
  passed_evts = totalGenFilterInfo_.sumPassWeights();
  filter_eff = totalGenFilterInfo_.filterEfficiency(-1);
  filter_eff_uncert = totalGenFilterInfo_.filterEfficiencyError(-1);
  
  outT->Fill();

}

// Define as a plugin
DEFINE_FWK_MODULE(GenFilterEfficiencySaver);
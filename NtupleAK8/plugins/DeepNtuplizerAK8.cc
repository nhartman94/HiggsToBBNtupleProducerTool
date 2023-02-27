/*
 * DeepNtuplizerAK8AK8.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DeepNTuples/NtupleCommons/interface/TreeWriter.h"

#include "DeepNTuples/NtupleAK8/interface/EvtInfoFiller.h"
#include "DeepNTuples/NtupleAK8/interface/JetInfoFillerAK8.h"
#include "DeepNTuples/NtupleAK8/interface/FatJetInfoFiller.h"
#include "DeepNTuples/NtupleAK8/interface/PFCandidateFiller.h"
#include "DeepNTuples/NtupleAK8/interface/TrackFiller.h"
#include "DeepNTuples/NtupleAK8/interface/SVFiller.h"


using namespace deepntuples;

class DeepNtuplizerAK8 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DeepNtuplizerAK8(const edm::ParameterSet&);
  ~DeepNtuplizerAK8();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;

  edm::Service<TFileService> fs;
  TreeWriter *treeWriter = nullptr;

  NtupleBase* addModule(NtupleBase *m){
    modules_.push_back(m);
    return m;
  }
  std::vector<NtupleBase*> modules_;
};

DeepNtuplizerAK8::DeepNtuplizerAK8(const edm::ParameterSet& iConfig) {

  // read configuration parameters
  const double jetR = iConfig.getParameter<double>("jetR");
  const double candidatePtMin = iConfig.getParameter<double>("candidatePtMin");

  std::cout << "jetR = " << jetR << std::endl;

  // register modules
  EvtInfoFiller *evtinfo = new EvtInfoFiller("", jetR);
  addModule(evtinfo);

  JetInfoFillerAK8 *jetinfo = new JetInfoFillerAK8("", jetR);
  addModule(jetinfo);

  // FatJetInfoFiller *fjinfo = new FatJetInfoFiller("", jetR);
  // addModule(fjinfo);

  // PFCandidateFiller *pfcands = new PFCandidateFiller("", jetR, candidatePtMin);
  // addModule(pfcands);

  TrackFiller *tracks = new TrackFiller("", jetR, candidatePtMin);
  addModule(tracks);

  // SVFiller *sv = new SVFiller("", jetR);
  // addModule(sv);

  // read config and init modules
  for(auto& m: modules_)
    m->readConfig(iConfig, consumesCollector());

}

DeepNtuplizerAK8::~DeepNtuplizerAK8()
{
  for(auto *m : modules_)
    delete m;
}


// ------------ method called for each event  ------------
void DeepNtuplizerAK8::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  for(auto *m : modules_){
    m->readEvent(iEvent, iSetup);
  }

  bool write_ = true;

  for (auto *m : modules_){
    if (!m->fillBranches()){
      std::cout << "Setting write_ to false\n";
      write_ = false;
      break;
    }
  }

  if (write_) {
    treeWriter->fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void DeepNtuplizerAK8::beginJob() {
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
        "TFile Service is not registered in cfg file" );
  }
  treeWriter = new TreeWriter(fs->make<TTree>("tree" ,"tree"));

  for(auto *m : modules_)
    m->initBranches(treeWriter);

}

// ------------ method called once each job just after ending the event loop  ------------
void DeepNtuplizerAK8::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DeepNtuplizerAK8::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizerAK8);

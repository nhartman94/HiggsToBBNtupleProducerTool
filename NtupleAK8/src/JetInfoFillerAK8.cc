/*
 * JetInfoAK8.cc
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleAK8/interface/JetInfoFillerAK8.h"

#include <algorithm>


namespace deepntuples {

void JetInfoFillerAK8::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) {

  jetToken_ = cc.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

  minPt_ = iConfig.getUntrackedParameter<double>("jetPtMin", 150);
  maxPt_ = iConfig.getUntrackedParameter<double>("jetPtMax", -1);
  maxAbsEta_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 2.4);
  btag_discriminators_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminators");

}

void JetInfoFillerAK8::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

//  iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);
//  iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);

  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
  flavorDef.setGenParticles(*genParticlesHandle);

//  iEvent.getByToken(muonsToken_, muonsHandle);
//  iEvent.getByToken(electronsToken_, electronsHandle);

  // This used to be in the steering script
  iEvent.getByToken(jetToken_, jets);

}

bool JetInfoFillerAK8::fill() {

  // This used to be in the steering script
  for (unsigned jetidx=0; jetidx<jets->size(); ++jetidx){
    
    const auto jet = jets->at(jetidx).correctedJet("Uncorrected"); // undo the JECs
    JetHelper jet_helper(&jet);

    // jet selection
    if ( (jet.pt() < minPt_) || (maxPt_ > 0 && jet.pt() > maxPt_)  ||  (std::abs(jet.eta()) > maxAbsEta_) )
      continue;

    // Now the code that was in the OG script... except with fillMulti instead of fill
    data.fillMulti<unsigned>("jet_no", jetidx);

    // truth labels
    float gen_pt = jet.genJet() ? jet.genJet()->pt() : 0;
    data.fillMulti<float>("gen_pt", gen_pt);
    data.fillMulti<float>("Delta_gen_pt", gen_pt - jet.correctedJet("Uncorrected").pt());

    auto flavor = flavorDef.jet_flavour(jet);
    data.fillMulti<int>("isB", flavor==JetFlavor::B);
    data.fillMulti<int>("isBB", flavor==JetFlavor::BB);
    data.fillMulti<int>("isLeptonicB", flavor==JetFlavor::LeptonicB);
    data.fillMulti<int>("isLeptonicB_C", flavor==JetFlavor::LeptonicB_C);
    data.fillMulti<int>("isC", flavor==JetFlavor::C);
    data.fillMulti<int>("isUD", flavor==JetFlavor::UD);
    data.fillMulti<int>("isS", flavor==JetFlavor::S);
    data.fillMulti<int>("isG", flavor==JetFlavor::G);
    data.fillMulti<int>("isUndefined", flavor==JetFlavor::UNDEFINED);

    // jet variables
    data.fillMulti<float>("jet_pt", jet.correctedJet("Uncorrected").pt());
    data.fillMulti<float>("jet_corr_pt", jet.pt());
    data.fillMulti<float>("jet_eta", jet.eta());
    data.fillMulti<float>("jet_phi", jet.phi());

    // jet id
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    bool jet_looseId_ = true;
    bool jet_tightId_ = true;
    try{
      float NHF  = jet.neutralHadronEnergyFraction();
      float NEMF = jet.neutralEmEnergyFraction();
      float CHF  = jet.chargedHadronEnergyFraction();
  //    float MUF  = jet.muonEnergyFraction();
      float CEMF = jet.chargedEmEnergyFraction();
      float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
      float NumNeutralParticles = jet.neutralMultiplicity();
      float CHM      = jet.chargedMultiplicity();

      jet_looseId_ = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet.eta())>2.4) && abs(jet.eta())<=2.7) ||
          (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(jet.eta())>2.7 && abs(jet.eta())<=3.0 ) ||
          (NEMF<0.90 && NumNeutralParticles>10 && abs(jet.eta())>3.0 );

      jet_tightId_ = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet.eta())>2.4) && abs(jet.eta())<=2.7 ) ||
          (NEMF<0.90 && NumNeutralParticles>2 && abs(jet.eta())>2.7 && abs(jet.eta())<=3.0) ||
          (NEMF<0.90 && NumNeutralParticles>10 && abs(jet.eta())>3.0);
    }catch(const cms::Exception &e){
      // energy fraction not supported on subjets/puppi?
    }

    data.fillMulti<float>("jet_looseId", jet_looseId_);
    data.fillMulti<float>("jet_tightId", jet_tightId_);

    for(const auto& disc : btag_discriminators_) {
      std::string name(disc);
      std::replace(name.begin(), name.end(), ':', '_');
      data.fillMulti<float>(name, catchInfs(jet.bDiscriminator(disc), -99));
    }

  }

  return true;
}

void JetInfoFillerAK8::book() {
  // event information
  data.addMulti<unsigned>("jet_no");

  // truth labels
  data.addMulti<float>("gen_pt");
  data.addMulti<float>("Delta_gen_pt");

  data.addMulti<int>("isB");
  data.addMulti<int>("isBB");
  data.addMulti<int>("isLeptonicB");
  data.addMulti<int>("isLeptonicB_C");
  data.addMulti<int>("isC");
  data.addMulti<int>("isUD");
  data.addMulti<int>("isS");
  data.addMulti<int>("isG");
  data.addMulti<int>("isUndefined");

  // jet variables
  data.addMulti<float>("jet_pt");
  data.addMulti<float>("jet_corr_pt");
  data.addMulti<float>("jet_eta");
  data.addMulti<float>("jet_phi");

  // jet id
  data.addMulti<float>("jet_looseId");
  data.addMulti<float>("jet_tightId");

  for(auto name : btag_discriminators_) {
    std::replace(name.begin(), name.end(), ':', '_');
    data.addMulti<float>(name);
  }
}

} /* namespace deepntuples */


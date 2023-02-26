/*
 * EvtInfoFiller.cc
 *
 *  Created on: Feb 24, 2023
 *      Author: nmh 
 */

#include "DeepNTuples/NtupleAK8/interface/EvtInfoFiller.h"

#include <algorithm>


using namespace deepntuples;
    
void EvtInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  puToken_ = cc.consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puInfo"));
  rhoToken_ = cc.consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"));
}


void EvtInfoFiller::readEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(puToken_, puInfo);
  iEvent.getByToken(rhoToken_, rhoInfo);
  event_ = iEvent.id().event();
}

bool EvtInfoFiller::fill() {
  
  // pv selection
  if (vertices->empty()) return false;

  // event information
  data.fill<unsigned>("event_no", event_);
  
  data.fill<float>("npv", vertices->size());
  data.fill<float>("rho", *rhoInfo);
  for (const auto &v : *puInfo) {
    int bx = v.getBunchCrossing();
    if (bx == 0) {
      data.fill<float>("ntrueInt", v.getTrueNumInteractions());
    }
  }

  return true;
}

void EvtInfoFiller::book() {

  // event information
  data.add<unsigned>("event_no", 0);
  data.add<float>("npv", 0);
  data.add<float>("rho", 0);
  data.add<float>("ntrueInt", 0);

}
    

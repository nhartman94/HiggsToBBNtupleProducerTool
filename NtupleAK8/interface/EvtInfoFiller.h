/*
 * EvtInfoFiller.h
 *
 *  Created on: Feb 24, 2023
 *      Author: nmh 
*/

#ifndef NTUPLEAK8_INTERFACE_EVTINFOFILLER_H_
#define NTUPLEAK8_INTERFACE_EVTINFOFILLER_H_

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class EvtInfoFiller: public NtupleBase {
public:
  EvtInfoFiller() : EvtInfoFiller("") {}
  EvtInfoFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR) {}
  virtual ~EvtInfoFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill() override;

private:
  
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
  edm::EDGetTokenT<double> rhoToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
  edm::Handle<double> rhoInfo;

  unsigned event_ = 0;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;

};

} /* namespace deepntuples */

#endif /* NTUPLEAK8_INTERFACE_EVTINFOFILLER_H_ */
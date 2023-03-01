/*
 * PFCandidateFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLEAK8_INTERFACE_PFCANDIDATEFILLER_H_
#define NTUPLEAK8_INTERFACE_PFCANDIDATEFILLER_H_

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class PFCandidateFiller: public NtupleBase {
public:
  PFCandidateFiller() : PFCandidateFiller("", 0.8, -1) {}
  PFCandidateFiller(std::string branchName, double jetR=0.8, double minPt=-1) : NtupleBase(branchName, jetR, minPt) {}
  virtual ~PFCandidateFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  // virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;
  virtual bool fill() override;

private:

  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  edm::Handle<edm::View<pat::Jet>> jets;
  
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;
};

} /* namespace deepntuples */

#endif /* NTUPLEAK8_INTERFACE_PFCANDIDATEFILLER_H_ */

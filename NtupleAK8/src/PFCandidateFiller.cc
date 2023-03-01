/*
 * PFCandidateFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleAK8/interface/PFCandidateFiller.h"

namespace deepntuples {

void PFCandidateFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  jetToken_ = cc.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCandidateFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(jetToken_, jets);
}

void PFCandidateFiller::book() {

  data.addMulti<int>("n_pfcands");
  data.addMulti<float>("npfcands");

  // basic kinematics, valid for both charged and neutral
  data.addMultiMulti<float>("pfcand_ptrel");
  data.addMultiMulti<float>("pfcand_erel");
  data.addMultiMulti<float>("pfcand_phirel");
  data.addMultiMulti<float>("pfcand_etarel");
  data.addMultiMulti<float>("pfcand_deltaR");
  data.addMultiMulti<float>("pfcand_puppiw");
  data.addMultiMulti<float>("pfcand_mass");

  data.addMultiMulti<float>("pfcand_drminsv");
  data.addMultiMulti<float>("pfcand_drsubjet1");
  data.addMultiMulti<float>("pfcand_drsubjet2");

  data.addMultiMulti<float>("pfcand_charge");
  data.addMultiMulti<float>("pfcand_isMu");
  data.addMultiMulti<float>("pfcand_isEl");
  data.addMultiMulti<float>("pfcand_isGamma");
  data.addMultiMulti<float>("pfcand_isChargedHad");
  data.addMultiMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMultiMulti<float>("pfcand_hcalFrac");

  // for charged
  data.addMultiMulti<float>("pfcand_VTX_ass");
  data.addMultiMulti<float>("pfcand_fromPV");
  data.addMultiMulti<float>("pfcand_lostInnerHits");

  // impact parameters
  data.addMultiMulti<float>("pfcand_dz");
  data.addMultiMulti<float>("pfcand_dzsig");
  data.addMultiMulti<float>("pfcand_dxy");
  data.addMultiMulti<float>("pfcand_dxysig");

}

bool PFCandidateFiller::fill() {

  for (unsigned jetidx=0; jetidx<jets->size(); ++jetidx){
    const auto jet = jets->at(jetidx).correctedJet("Uncorrected"); // undo the JECs
    JetHelper jet_helper(&jet);
    const auto& jetConstituents = jet_helper.getJetConstituents();
    int nConstituents = 0;
    float etasign = jet.eta()>0 ? 1 : -1;

    for (const auto *pfcand : jetConstituents){

      if (pfcand->pt() < minPt_) continue;
      nConstituents += 1;
      // basic kinematics, valid for both charged and neutral
      data.fillMultiMulti<float>("pfcand_ptrel", pfcand->pt()/jet.pt(), jetidx);
      data.fillMultiMulti<float>("pfcand_erel", pfcand->energy()/jet.energy(), jetidx);
      data.fillMultiMulti<float>("pfcand_phirel", reco::deltaPhi(*pfcand, jet), jetidx);
      data.fillMultiMulti<float>("pfcand_etarel", etasign * (pfcand->eta() - jet.eta()), jetidx);
      data.fillMultiMulti<float>("pfcand_deltaR", reco::deltaR(*pfcand, jet), jetidx);
      data.fillMultiMulti<float>("pfcand_puppiw", pfcand->puppiWeight(), jetidx);
      data.fillMultiMulti<float>("pfcand_mass", pfcand->mass(), jetidx);
      double minDR = 999;
      for (const auto &sv : *SVs){
        double dr = reco::deltaR(*pfcand, sv);
        if (dr < minDR) minDR = dr;
      }
      data.fillMultiMulti<float>("pfcand_drminsv", minDR==999 ? -1 : minDR, jetidx);
      const auto& subjets = jet_helper.getSubJets();
      data.fillMultiMulti<float>("pfcand_drsubjet1", subjets.size()>0 ? reco::deltaR(*pfcand, *subjets.at(0)) : -1, jetidx);
      data.fillMultiMulti<float>("pfcand_drsubjet2", subjets.size()>1 ? reco::deltaR(*pfcand, *subjets.at(1)) : -1, jetidx);

      data.fillMultiMulti<float>("pfcand_charge", pfcand->charge(), jetidx);
      data.fillMultiMulti<float>("pfcand_isEl", std::abs(pfcand->pdgId())==11, jetidx);
      data.fillMultiMulti<float>("pfcand_isMu", std::abs(pfcand->pdgId())==13, jetidx);
      data.fillMultiMulti<float>("pfcand_isGamma", std::abs(pfcand->pdgId())==22, jetidx);
      data.fillMultiMulti<float>("pfcand_isChargedHad", std::abs(pfcand->pdgId())==211, jetidx);
      data.fillMultiMulti<float>("pfcand_isNeutralHad", std::abs(pfcand->pdgId())==130, jetidx);

      // for neutral
      data.fillMultiMulti<float>("pfcand_hcalFrac", pfcand->hcalFraction(), jetidx);

      // for charged
      data.fillMultiMulti<float>("pfcand_VTX_ass", pfcand->pvAssociationQuality(), jetidx);
      data.fillMultiMulti<float>("pfcand_fromPV", pfcand->fromPV(), jetidx);
      data.fillMultiMulti<float>("pfcand_lostInnerHits", pfcand->lostInnerHits(), jetidx);

      // impact parameters
      data.fillMultiMulti<float>("pfcand_dz", catchInfs(pfcand->dz()), jetidx);
      data.fillMultiMulti<float>("pfcand_dxy", catchInfs(pfcand->dxy()), jetidx);
      //if( pfcand->hasTrackDetails() ) { // 101X
      if( true ) { // 80X
        data.fillMultiMulti<float>("pfcand_dzsig", catchInfsAndBound(pfcand->dz()/pfcand->dzError(),0,-2000,2000), jetidx);
        data.fillMultiMulti<float>("pfcand_dxysig", catchInfsAndBound(pfcand->dxy()/pfcand->dxyError(),0,-2000,2000), jetidx);
      }
      else {
        data.fillMultiMulti<float>("pfcand_dzsig", -1, jetidx);
        data.fillMultiMulti<float>("pfcand_dxysig", -1, jetidx);
      }

    }
    data.fillMulti<int>("n_pfcands", nConstituents);
    data.fillMulti<float>("npfcands", nConstituents);

  }
  return true;
}

} /* namespace deepntuples */

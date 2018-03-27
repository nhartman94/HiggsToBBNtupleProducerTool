/*
 * TrackFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/NtupleAK8/interface/TrackFiller.h"

#include "DeepNTuples/NtupleCommons/interface/sorting_modules.h"
#include "DeepNTuples/NtupleCommons/interface/InfinityCatcher.h"
namespace deepntuples {

void TrackFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void TrackFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void TrackFiller::book() {

  data.add<int>("n_tracks", 0);
  data.add<float>("ntracks", 0);

  // basic kinematics
  data.addMulti<float>("track_ptrel");
  data.addMulti<float>("track_erel");
  data.addMulti<float>("track_phirel");
  data.addMulti<float>("track_etarel");
  data.addMulti<float>("track_deltaR");
  data.addMulti<float>("track_puppiw");
  data.addMulti<float>("track_pt");
  data.addMulti<float>("track_mass");

  data.addMulti<float>("track_drminsv");
  data.addMulti<float>("track_drsubjet1");
  data.addMulti<float>("track_drsubjet2");

  data.addMulti<float>("track_charge");
  data.addMulti<float>("track_isMu");
  data.addMulti<float>("track_isEl");
  data.addMulti<float>("track_isChargedHad");

  // for charged
  data.addMulti<float>("track_VTX_ass");
  data.addMulti<float>("track_fromPV");
  data.addMulti<float>("track_lostInnerHits");

  // impact parameters
  data.addMulti<float>("track_dz");
  data.addMulti<float>("track_dzsig");
  data.addMulti<float>("track_dxy");
  data.addMulti<float>("track_dxysig");

  // track quality
  data.addMulti<float>("track_normchi2");
  data.addMulti<float>("track_quality");

  // track covariance
  data.addMulti<float>("track_dptdpt");
  data.addMulti<float>("track_detadeta");
  data.addMulti<float>("track_dphidphi");
  data.addMulti<float>("track_dxydxy");
  data.addMulti<float>("track_dzdz");
  data.addMulti<float>("track_dxydz");
  data.addMulti<float>("track_dphidxy");
  data.addMulti<float>("track_dlambdadz");

  // track btag info
  data.addMulti<float>("trackBTag_Momentum");
  data.addMulti<float>("trackBTag_Eta");
  data.addMulti<float>("trackBTag_EtaRel");
  data.addMulti<float>("trackBTag_PtRel");
  data.addMulti<float>("trackBTag_PPar");
  data.addMulti<float>("trackBTag_DeltaR");
  data.addMulti<float>("trackBTag_PtRatio");
  data.addMulti<float>("trackBTag_PParRatio");
  data.addMulti<float>("trackBTag_Sip2dVal");
  data.addMulti<float>("trackBTag_Sip2dSig");
  data.addMulti<float>("trackBTag_Sip3dVal");
  data.addMulti<float>("trackBTag_Sip3dSig");
  data.addMulti<float>("trackBTag_JetDistVal");
//  data.addMulti<float>("trackBTag_JetDistSig"); // always gives 0


}

bool TrackFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const pat::PackedCandidate*> chargedPFCands;
  std::unordered_map<const pat::PackedCandidate*, TrackInfoBuilder> trackInfoMap;
  std::unordered_map<const pat::PackedCandidate*, double> drMinSvMap;
  std::unordered_map<const pat::PackedCandidate*, double> ptRelMap;

  std::vector<sorting::sortingClass<size_t>> sortedcharged;

  const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();

  unsigned int i = 0;
  for (const auto * pfcand : jet_helper.getJetConstituents()){
    if (!pfcand) continue;
    if (pfcand->pt() < minPt_) continue;
    if (pfcand->charge() != 0) {
      chargedPFCands.push_back(pfcand);
      trackInfoMap[pfcand];
      trackInfoMap[pfcand].buildTrackInfo(builder_, *pfcand, jet, vertices->at(0));
      drMinSvMap[pfcand];
      double minDR = 0.8;
      for (const auto &sv : *SVs){
	double dr = reco::deltaR(*pfcand, sv);
	if (dr < minDR) minDR = dr;
      }
      drMinSvMap[pfcand] = minDR;

      sortedcharged.push_back(sorting::sortingClass<size_t>
			      (i, trackInfoMap.at(pfcand).getTrackSip2dSigRaw(),
			       -minDR, pfcand->pt()/jet_uncorr_pt));
      i++;
    }
  }

  // sort by ABCInv
  std::sort(sortedcharged.begin(),sortedcharged.end(),sorting::sortingClass<size_t>::compareByABCInv);
  std::vector<size_t> sortedchargedindices;
  sortedchargedindices=sorting::invertSortingVector(sortedcharged);

  // sort by Sip2d significance
  //std::sort(chargedPFCands.begin(), chargedPFCands.end(), [&](const pat::PackedCandidate *p1, const pat::PackedCandidate *p2){
  //  return trackInfoMap.at(p1).getTrackSip2dSig() > trackInfoMap.at(p2).getTrackSip2dSig();
  //});

  data.fill<int>("n_tracks", chargedPFCands.size());
  data.fill<float>("ntracks", chargedPFCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  //  for (const auto *cpf : chargedPFCands){
  for (i = 0; i < chargedPFCands.size(); i++){
    //std::cout << "sortedchargedindices.at(" << i << ") = " << sortedchargedindices.at(i) << std::endl;
    //std::cout << "sortedcharged.at(" << i << ") = " << sortedcharged.at(i).get() << std::endl;

    //const auto *cpf = chargedPFCands.at(sortedchargedindices.at(i));
    const auto *cpf = chargedPFCands.at(sortedcharged.at(i).get());

    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("track_ptrel", cpf->pt()/jet.pt());
    data.fillMulti<float>("track_erel", cpf->energy()/jet.energy());
    data.fillMulti<float>("track_phirel", reco::deltaPhi(*cpf, jet));
    data.fillMulti<float>("track_etarel", etasign * (cpf->eta() - jet.eta()));
    data.fillMulti<float>("track_deltaR", reco::deltaR(*cpf, jet));
    data.fillMulti<float>("track_puppiw", cpf->puppiWeight());
    data.fillMulti<float>("track_pt", cpf->pt());
    data.fillMulti<float>("track_mass", cpf->mass());

    data.fillMulti<float>("track_drminsv", catchInfsAndBound(drMinSvMap.at(cpf),0,-0.8,0,-0.8));

    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("track_drsubjet1", subjets.size()>0 ? reco::deltaR(*cpf, *subjets.at(0)) : -1);
    data.fillMulti<float>("track_drsubjet2", subjets.size()>1 ? reco::deltaR(*cpf, *subjets.at(1)) : -1);

    data.fillMulti<float>("track_charge", cpf->charge());
    data.fillMulti<float>("track_isEl", std::abs(cpf->pdgId())==11);
    data.fillMulti<float>("track_isMu", std::abs(cpf->pdgId())==13);
    data.fillMulti<float>("track_isChargedHad", std::abs(cpf->pdgId())==211);

    // for charged
    data.fillMulti<float>("track_VTX_ass", cpf->pvAssociationQuality());
    data.fillMulti<float>("track_fromPV", cpf->fromPV());
    data.fillMulti<float>("track_lostInnerHits", cpf->lostInnerHits());

    // impact parameters
    data.fillMulti<float>("track_dxy", catchInfs(cpf->dxy()));
    data.fillMulti<float>("track_dz", catchInfs(cpf->dz()));

    //if( cpf->hasTrackDetails() ) { // 101X
    if( true ) { // 80X
      const auto &trk = cpf->pseudoTrack();
      data.fillMulti<float>("track_dxysig", catchInfsAndBound(cpf->dxy()/cpf->dxyError(),0,-2000,2000));
      data.fillMulti<float>("track_dzsig", catchInfsAndBound(cpf->dz()/cpf->dzError(),0,-2000,2000));   
      data.fillMulti<float>("track_normchi2", catchInfsAndBound(trk.normalizedChi2(),300,-1,300));
      data.fillMulti<float>("track_quality", trk.qualityMask());
      // track covariance
      auto cov = [&](unsigned i, unsigned j) {
	return catchInfs(trk.covariance(i, j));
      };
      data.fillMulti<float>("track_dptdpt", cov(0,0));
      data.fillMulti<float>("track_detadeta", cov(1,1));
      data.fillMulti<float>("track_dphidphi", cov(2,2));
      data.fillMulti<float>("track_dxydxy", cov(3,3));
      data.fillMulti<float>("track_dzdz", cov(4,4));
      data.fillMulti<float>("track_dxydz", cov(3,4));
      data.fillMulti<float>("track_dphidxy", cov(2,3));
      data.fillMulti<float>("track_dlambdadz", cov(1,4));
    }
    else {
      data.fillMulti<float>("track_dxysig", -1);
      data.fillMulti<float>("track_dzsig", -1);
      data.fillMulti<float>("track_normchi2", -1);
      data.fillMulti<float>("track_quality",(1 << reco::TrackBase::loose));
      data.fillMulti<float>("track_dptdpt", 0);
      data.fillMulti<float>("track_detadeta", 0);
      data.fillMulti<float>("track_dphidphi", 0);
      data.fillMulti<float>("track_dxydxy", 0);
      data.fillMulti<float>("track_dzdz", 0);
      data.fillMulti<float>("track_dxydz", 0);
      data.fillMulti<float>("track_dphidxy", 0);
      data.fillMulti<float>("track_dlambdadz", 0);
    }

    const auto &trkinfo = trackInfoMap.at(cpf);
    data.fillMulti<float>("trackBTag_Momentum", trkinfo.getTrackMomentum());
    data.fillMulti<float>("trackBTag_Eta", trkinfo.getTrackEta());
    data.fillMulti<float>("trackBTag_EtaRel", trkinfo.getTrackEtaRel());
    data.fillMulti<float>("trackBTag_PtRel", trkinfo.getTrackPtRel());
    data.fillMulti<float>("trackBTag_PPar", trkinfo.getTrackPPar());
    data.fillMulti<float>("trackBTag_DeltaR", trkinfo.getTrackDeltaR());
    data.fillMulti<float>("trackBTag_PtRatio", trkinfo.getTrackPtRatio());
    data.fillMulti<float>("trackBTag_PParRatio", trkinfo.getTrackPParRatio());
    data.fillMulti<float>("trackBTag_Sip2dVal", trkinfo.getTrackSip2dVal());
    data.fillMulti<float>("trackBTag_Sip2dSig", trkinfo.getTrackSip2dSig());
    data.fillMulti<float>("trackBTag_Sip3dVal", trkinfo.getTrackSip3dVal());
    data.fillMulti<float>("trackBTag_Sip3dSig", trkinfo.getTrackSip3dSig());
    data.fillMulti<float>("trackBTag_JetDistVal", trkinfo.getTrackJetDistVal());
//    data.fillMulti<float>("trackBTag_JetDistSig", trkinfo.getTrackJetDistSig());

  }


  return true;
}

} /* namespace deepntuples */

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

  jetToken_ = cc.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));

  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void TrackFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  iEvent.getByToken(jetToken_, jets);  

  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void TrackFiller::book() {

  data.addMulti<int>("n_tracks");
  data.addMulti<float>("ntracks");

  // basic kinematics
  data.addMultiMulti<float>("track_ptrel");
  data.addMultiMulti<float>("track_erel");
  data.addMultiMulti<float>("track_phirel");
  data.addMultiMulti<float>("track_etarel");
  data.addMultiMulti<float>("track_deltaR");
  data.addMultiMulti<float>("track_puppiw");
  data.addMultiMulti<float>("track_pt");
  data.addMultiMulti<float>("track_mass");

  data.addMultiMulti<float>("track_drminsv");
  data.addMultiMulti<float>("track_drsubjet1");
  data.addMultiMulti<float>("track_drsubjet2");

  data.addMultiMulti<float>("track_charge");
  data.addMultiMulti<float>("track_isMu");
  data.addMultiMulti<float>("track_isEl");
  data.addMultiMulti<float>("track_isChargedHad");

  // for charged
  data.addMultiMulti<float>("track_VTX_ass");
  data.addMultiMulti<float>("track_fromPV");
  data.addMultiMulti<float>("track_lostInnerHits");

  // impact parameters
  data.addMultiMulti<float>("track_dz");
  data.addMultiMulti<float>("track_dzsig");
  data.addMultiMulti<float>("track_dxy");
  data.addMultiMulti<float>("track_dxysig");

  // track quality
  data.addMultiMulti<float>("track_normchi2");
  data.addMultiMulti<float>("track_quality");

  // track covariance
  data.addMultiMulti<float>("track_dptdpt");
  data.addMultiMulti<float>("track_detadeta");
  data.addMultiMulti<float>("track_dphidphi");
  data.addMultiMulti<float>("track_dxydxy");
  data.addMultiMulti<float>("track_dzdz");
  data.addMultiMulti<float>("track_dxydz");
  data.addMultiMulti<float>("track_dphidxy");
  data.addMultiMulti<float>("track_dlambdadz");

  // track btag info
  data.addMultiMulti<float>("trackBTag_Momentum");
  data.addMultiMulti<float>("trackBTag_Eta");
  data.addMultiMulti<float>("trackBTag_EtaRel");
  data.addMultiMulti<float>("trackBTag_PtRel");
  data.addMultiMulti<float>("trackBTag_PPar");
  data.addMultiMulti<float>("trackBTag_DeltaR");
  data.addMultiMulti<float>("trackBTag_PtRatio");
  data.addMultiMulti<float>("trackBTag_PParRatio");
  data.addMultiMulti<float>("trackBTag_Sip2dVal");
  data.addMultiMulti<float>("trackBTag_Sip2dSig");
  data.addMultiMulti<float>("trackBTag_Sip3dVal");
  data.addMultiMulti<float>("trackBTag_Sip3dSig");
  data.addMultiMulti<float>("trackBTag_JetDistVal");
//  data.addMulti<float>("trackBTag_JetDistSig"); // always gives 0

}

bool TrackFiller::fill() {

  std::vector<const pat::PackedCandidate*> chargedPFCands;
  std::unordered_map<const pat::PackedCandidate*, TrackInfoBuilder> trackInfoMap;
  std::unordered_map<const pat::PackedCandidate*, double> drMinSvMap;
  std::unordered_map<const pat::PackedCandidate*, double> ptRelMap;

  std::vector<sorting::sortingClass<size_t>> sortedcharged;

  for (unsigned jetidx=0; jetidx<jets->size(); ++jetidx){
    
    const auto jet = jets->at(jetidx).correctedJet("Uncorrected"); // undo the JECs
    JetHelper jet_helper(&jet);

    const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();
    chargedPFCands = {};
    sortedcharged = {};
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

    data.fillMulti<int>("n_tracks", chargedPFCands.size());
    data.fillMulti<float>("ntracks", chargedPFCands.size());

    float etasign = jet.eta()>0 ? 1 : -1;

    //  for (const auto *cpf : chargedPFCands){
    for (i = 0; i < chargedPFCands.size(); i++){
      //std::cout << "sortedchargedindices.at(" << i << ") = " << sortedchargedindices.at(i) << std::endl;
      //std::cout << "sortedcharged.at(" << i << ") = " << sortedcharged.at(i).get() << std::endl;

      //const auto *cpf = chargedPFCands.at(sortedchargedindices.at(i));
      const auto *cpf = chargedPFCands.at(sortedcharged.at(i).get());

      // basic kinematics, valid for both charged and neutral
      data.fillMultiMulti<float>("track_ptrel", cpf->pt()/jet.pt(), jetidx);
      data.fillMultiMulti<float>("track_erel", cpf->energy()/jet.energy(), jetidx);
      data.fillMultiMulti<float>("track_phirel", reco::deltaPhi(*cpf, jet), jetidx);
      data.fillMultiMulti<float>("track_etarel", etasign * (cpf->eta() - jet.eta()), jetidx);
      data.fillMultiMulti<float>("track_deltaR", reco::deltaR(*cpf, jet), jetidx);
      data.fillMultiMulti<float>("track_puppiw", cpf->puppiWeight(), jetidx);
      data.fillMultiMulti<float>("track_pt", cpf->pt(), jetidx);
      data.fillMultiMulti<float>("track_mass", cpf->mass(), jetidx);

      data.fillMultiMulti<float>("track_drminsv", catchInfsAndBound(drMinSvMap.at(cpf),0,-0.8,0,-0.8), jetidx);

      const auto& subjets = jet_helper.getSubJets();
      data.fillMultiMulti<float>("track_drsubjet1", subjets.size()>0 ? reco::deltaR(*cpf, *subjets.at(0)) : -1, jetidx);
      data.fillMultiMulti<float>("track_drsubjet2", subjets.size()>1 ? reco::deltaR(*cpf, *subjets.at(1)) : -1, jetidx);

      data.fillMultiMulti<float>("track_charge", cpf->charge(), jetidx);
      data.fillMultiMulti<float>("track_isEl", std::abs(cpf->pdgId())==11, jetidx);
      data.fillMultiMulti<float>("track_isMu", std::abs(cpf->pdgId())==13, jetidx);
      data.fillMultiMulti<float>("track_isChargedHad", std::abs(cpf->pdgId())==211, jetidx);

      // for charged
      data.fillMultiMulti<float>("track_VTX_ass", cpf->pvAssociationQuality(), jetidx);
      data.fillMultiMulti<float>("track_fromPV", cpf->fromPV(), jetidx);
      data.fillMultiMulti<float>("track_lostInnerHits", cpf->lostInnerHits(), jetidx);

      // impact parameters
      data.fillMultiMulti<float>("track_dxy", catchInfs(cpf->dxy()), jetidx);
      data.fillMultiMulti<float>("track_dz", catchInfs(cpf->dz()), jetidx);

      //if( cpf->hasTrackDetails() ) { // 101X
      if( true ) { // 80X
        const auto &trk = cpf->pseudoTrack();
        data.fillMultiMulti<float>("track_dxysig", catchInfsAndBound(cpf->dxy()/cpf->dxyError(),0,-2000,2000), jetidx);
        data.fillMultiMulti<float>("track_dzsig", catchInfsAndBound(cpf->dz()/cpf->dzError(),0,-2000,2000), jetidx);   
        data.fillMultiMulti<float>("track_normchi2", catchInfsAndBound(trk.normalizedChi2(),300,-1,300), jetidx);
        data.fillMultiMulti<float>("track_quality", trk.qualityMask(), jetidx);
        // track covariance
        auto cov = [&](unsigned i, unsigned j) {
    return catchInfs(trk.covariance(i, j));
        };
        data.fillMultiMulti<float>("track_dptdpt", cov(0,0), jetidx);
        data.fillMultiMulti<float>("track_detadeta", cov(1,1), jetidx);
        data.fillMultiMulti<float>("track_dphidphi", cov(2,2), jetidx);
        data.fillMultiMulti<float>("track_dxydxy", cov(3,3), jetidx);
        data.fillMultiMulti<float>("track_dzdz", cov(4,4), jetidx);
        data.fillMultiMulti<float>("track_dxydz", cov(3,4), jetidx);
        data.fillMultiMulti<float>("track_dphidxy", cov(2,3), jetidx);
        data.fillMultiMulti<float>("track_dlambdadz", cov(1,4), jetidx);
      }
      else {
        data.fillMultiMulti<float>("track_dxysig", -1, jetidx);
        data.fillMultiMulti<float>("track_dzsig", -1, jetidx);
        data.fillMultiMulti<float>("track_normchi2", -1, jetidx);
        data.fillMultiMulti<float>("track_quality",(1 << reco::TrackBase::loose), jetidx);
        data.fillMultiMulti<float>("track_dptdpt", 0, jetidx);
        data.fillMultiMulti<float>("track_detadeta", 0, jetidx);
        data.fillMultiMulti<float>("track_dphidphi", 0, jetidx);
        data.fillMultiMulti<float>("track_dxydxy", 0, jetidx);
        data.fillMultiMulti<float>("track_dzdz", 0, jetidx);
        data.fillMultiMulti<float>("track_dxydz", 0, jetidx);
        data.fillMultiMulti<float>("track_dphidxy", 0, jetidx);
        data.fillMultiMulti<float>("track_dlambdadz", 0, jetidx);
      }

      const auto &trkinfo = trackInfoMap.at(cpf);
      data.fillMultiMulti<float>("trackBTag_Momentum", trkinfo.getTrackMomentum(), jetidx);
      data.fillMultiMulti<float>("trackBTag_Eta", trkinfo.getTrackEta(), jetidx);
      data.fillMultiMulti<float>("trackBTag_EtaRel", trkinfo.getTrackEtaRel(), jetidx);
      data.fillMultiMulti<float>("trackBTag_PtRel", trkinfo.getTrackPtRel(), jetidx);
      data.fillMultiMulti<float>("trackBTag_PPar", trkinfo.getTrackPPar(), jetidx);
      data.fillMultiMulti<float>("trackBTag_DeltaR", trkinfo.getTrackDeltaR(), jetidx);
      data.fillMultiMulti<float>("trackBTag_PtRatio", trkinfo.getTrackPtRatio(), jetidx);
      data.fillMultiMulti<float>("trackBTag_PParRatio", trkinfo.getTrackPParRatio(), jetidx);
      data.fillMultiMulti<float>("trackBTag_Sip2dVal", trkinfo.getTrackSip2dVal(), jetidx);
      data.fillMultiMulti<float>("trackBTag_Sip2dSig", trkinfo.getTrackSip2dSig(), jetidx);
      data.fillMultiMulti<float>("trackBTag_Sip3dVal", trkinfo.getTrackSip3dVal(), jetidx);
      data.fillMultiMulti<float>("trackBTag_Sip3dSig", trkinfo.getTrackSip3dSig(), jetidx);
      data.fillMultiMulti<float>("trackBTag_JetDistVal", trkinfo.getTrackJetDistVal(), jetidx);
//    data.fillMulti<float>("trackBTag_JetDistSig", trkinfo.getTrackJetDistSig());

  }

  }


  return true;
}

} /* namespace deepntuples */

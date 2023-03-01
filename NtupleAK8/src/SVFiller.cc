/*
 * SVFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DeepNTuples/NtupleAK8/interface/SVFiller.h"

namespace deepntuples {

void SVFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  
  jetToken_ = cc.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));

  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void SVFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  iEvent.getByToken(jetToken_, jets);

  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
}

void SVFiller::book() {

  data.addMulti<int>("n_sv");
  data.addMulti<float>("nsv");

  // basic kinematics
  data.addMultiMulti<float>("sv_ptrel");
  data.addMultiMulti<float>("sv_erel");
  data.addMultiMulti<float>("sv_phirel");
  data.addMultiMulti<float>("sv_etarel");
  data.addMultiMulti<float>("sv_deltaR");
  data.addMultiMulti<float>("sv_pt");
  data.addMultiMulti<float>("sv_mass");

  // sv properties
  data.addMultiMulti<float>("sv_ntracks");
  data.addMultiMulti<float>("sv_chi2");
  data.addMultiMulti<float>("sv_ndf");
  data.addMultiMulti<float>("sv_normchi2");
  data.addMultiMulti<float>("sv_dxy");
  data.addMultiMulti<float>("sv_dxyerr");
  data.addMultiMulti<float>("sv_dxysig");
  data.addMultiMulti<float>("sv_d3d");
  data.addMultiMulti<float>("sv_d3derr");
  data.addMultiMulti<float>("sv_d3dsig");
  data.addMultiMulti<float>("sv_costhetasvpv");

}

bool SVFiller::fill() {

  std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;

  for (unsigned jetidx=0; jetidx<jets->size(); ++jetidx){
    
    const auto jet = jets->at(jetidx).correctedJet("Uncorrected"); // undo the JECs
    JetHelper jet_helper(&jet);

    jetSVs = {};
    for (const auto &sv : *SVs){
      if (reco::deltaR(sv, jet) < jetR_) {
        jetSVs.push_back(&sv);
      }
    }

    // sort by dxy significance
    const auto &pv = vertices->at(0);
    std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sv1, const reco::VertexCompositePtrCandidate *sv2){
      return vertexDxy(*sv1, pv).significance() > vertexDxy(*sv2, pv).significance();
    });


    data.fillMulti<int>("n_sv", jetSVs.size());
    data.fillMulti<float>("nsv", jetSVs.size());

    float etasign = jet.eta()>0 ? 1 : -1;

    for (const auto *sv : jetSVs){
      // basic kinematics
      data.fillMultiMulti<float>("sv_ptrel", sv->pt() / jet.pt(), jetidx);
      data.fillMultiMulti<float>("sv_erel", sv->energy() / jet.energy(), jetidx);
      data.fillMultiMulti<float>("sv_phirel", reco::deltaPhi(*sv, jet), jetidx);
      data.fillMultiMulti<float>("sv_etarel", etasign * (sv->eta() - jet.eta()), jetidx);
      data.fillMultiMulti<float>("sv_deltaR", catchInfsAndBound(std::fabs(reco::deltaR(*sv,jet))-0.5,0,-2,0), jetidx);
      data.fillMultiMulti<float>("sv_pt", sv->pt(), jetidx);
      data.fillMultiMulti<float>("sv_mass", sv->mass(), jetidx);

      // sv properties
      data.fillMultiMulti<float>("sv_ntracks", sv->numberOfDaughters(), jetidx);
      data.fillMultiMulti<float>("sv_chi2", sv->vertexChi2(), jetidx);
      data.fillMultiMulti<float>("sv_ndf", sv->vertexNdof(), jetidx);
      data.fillMultiMulti<float>("sv_normchi2", catchInfsAndBound(sv->vertexNormalizedChi2(), 1000, -1000, 1000), jetidx);

      const auto &dxy = vertexDxy(*sv, pv);
      data.fillMultiMulti<float>("sv_dxy", dxy.value(), jetidx);
      data.fillMultiMulti<float>("sv_dxyerr", dxy.error(), jetidx);
      data.fillMultiMulti<float>("sv_dxysig", catchInfsAndBound(dxy.significance(), 0,-1,800), jetidx);

      const auto &d3d = vertexD3d(*sv, pv);
      data.fillMultiMulti<float>("sv_d3d", d3d.value(), jetidx);
      data.fillMultiMulti<float>("sv_d3derr", d3d.error(), jetidx);
      data.fillMultiMulti<float>("sv_d3dsig", catchInfsAndBound(d3d.significance(), 0,-1,800), jetidx);
      data.fillMultiMulti<float>("sv_costhetasvpv", vertexDdotP(*sv, pv), jetidx);
    }
  }
  return true;
}


Measurement1D SVFiller::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D SVFiller::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float SVFiller::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}



} /* namespace deepntuples */

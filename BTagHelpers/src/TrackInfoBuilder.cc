#include <cmath>
#include "TVector3.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "DeepNTuples/NtupleCommons/interface/InfinityCatcher.h"
#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"

namespace deepntuples {

void TrackInfoBuilder::buildTrackInfo(const edm::ESHandle<TransientTrackBuilder> builder,
    const pat::PackedCandidate &pfcand, const pat::Jet &jet, const reco::Vertex& pv) {

  // DataFormats/BTauReco/interface/IPTagInfo.h

  math::XYZVector jetDir = jet.momentum().Unit();
  TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());

  GlobalVector jetXYZVector(jet.px(),jet.py(),jet.pz());
  //std::cout << "pv = " << pv.x() << pv.y() << pv.z() << std::endl;

    
  if(!pfcand.hasTrackDetails()) {
    TVector3 trackMom3(
		       pfcand.momentum().x(),
		       pfcand.momentum().y(),
		       pfcand.momentum().z()
		       );
    trackMomentum_=catchInfsAndBound(pfcand.p(),0,0 ,1000);
    trackEta_= catchInfsAndBound(pfcand.eta(),0,-5,5);
    trackEtaRel_=catchInfsAndBound(reco::btau::etaRel(jetDir, pfcand.momentum()),0,-5,15);
    trackPtRel_=catchInfsAndBound(trackMom3.Perp(jetDir3),0,-1,4);
    trackPPar_=catchInfsAndBound(jetDir.Dot(pfcand.momentum()),0,-1e5,1e5);
    trackDeltaR_=catchInfsAndBound(reco::deltaR(pfcand.momentum(), jetDir),0,-5,5);
    trackPtRatio_=catchInfsAndBound(trackMom3.Perp(jetDir3) / pfcand.p(),0,-1,4);
    trackPParRatio_=catchInfsAndBound(jetDir.Dot(pfcand.momentum()) / pfcand.p(),0,-10,100);
    trackSip2dVal_=0.;
    trackSip2dSig_=0.;
    trackSip3dVal_=0.;
    trackSip3dSig_=0.;
    trackJetDistVal_=0.;
    trackJetDistSig_=0.;
    return;
  }
  const auto &trk = pfcand.pseudoTrack();
  reco::TransientTrack transientTrack(builder->build(trk));
  Measurement1D meas_ip2d = IPTools::signedTransverseImpactParameter(transientTrack, jetXYZVector, pv).second;
  Measurement1D meas_ip3d = IPTools::signedImpactParameter3D(transientTrack, jetXYZVector, pv).second;
  Measurement1D jetdist = IPTools::jetTrackDistance(transientTrack, jetXYZVector, pv).second;
  math::XYZVector trackMom = trk.momentum();
  double trackMag = std::sqrt(trackMom.Mag2());

  TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());

  trackMomentum_ = catchInfsAndBound(trackMag,0,0 ,1000);
  trackEta_ = catchInfsAndBound(trackMom.Eta(),0,-5,5);
  trackEtaRel_ = catchInfsAndBound(reco::btau::etaRel(jetDir, trackMom),0,-5,15);
  trackPtRel_ = catchInfsAndBound(trackMom3.Perp(jetDir3),0,-1,4);
  trackPPar_ = catchInfsAndBound(jetDir.Dot(trackMom),0,-1e5,1e5);
  trackDeltaR_ = catchInfsAndBound(reco::deltaR(trackMom, jetDir),0,-5,5);
  trackPtRatio_ = catchInfsAndBound(trackMom3.Perp(jetDir3) / trackMag, 0, -1, 10);
  trackPParRatio_ = catchInfsAndBound(jetDir.Dot(trackMom) / trackMag,0,-10,100);
  trackSip2dVal_ = catchInfsAndBound(meas_ip2d.value(),0,-1,70);
  trackSip2dSig_ = catchInfsAndBound(meas_ip2d.significance(),0,-1,4e4);
  trackSip2dSigRaw_ = meas_ip2d.significance();
  trackSip3dVal_ = catchInfsAndBound(meas_ip3d.value(),0,-1,1e5);
  trackSip3dSig_ = catchInfsAndBound(meas_ip3d.significance(),0,-1,4e4);
  trackJetDistVal_ = catchInfsAndBound(jetdist.value(),0,-20,1);
  trackJetDistSig_ = catchInfsAndBound(jetdist.significance(),0,-1,1e5);
}

} /* namespace deepntuples */


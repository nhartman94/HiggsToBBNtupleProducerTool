/*
 * FatJetInfoFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleAK8/interface/FatJetInfoFiller.h"

namespace deepntuples {

void FatJetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  jetToken_ = cc.consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  fjTagInfoName = iConfig.getParameter<std::string>("fjTagInfoName");
  for (const auto &flv : iConfig.getUntrackedParameter<std::vector<unsigned>>("fjKeepFlavors", {})){
    keepFlavors_.push_back(static_cast<FatJetMatching::FatJetFlavor>(flv));
  }
  isQCDSample_ = iConfig.getUntrackedParameter<bool>("isQCDSample", false);
}

void FatJetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
  iEvent.getByToken(jetToken_, jets);
}

void FatJetInfoFiller::book() {
  // truth labels
  data.addMulti<int>("fj_label");
  data.addMulti<int>("fj_isTop");
  data.addMulti<int>("fj_isW");
  data.addMulti<int>("fj_isZ");
  data.addMulti<int>("fj_isH");
  data.addMulti<int>("fj_isQCD");

  data.addMulti<int>("label_Top_bcq");
  data.addMulti<int>("label_Top_bqq");
  data.addMulti<int>("label_Top_bc");
  data.addMulti<int>("label_Top_bq");

  data.addMulti<int>("label_W_cq");
  data.addMulti<int>("label_W_qq");

  data.addMulti<int>("label_Z_bb");
  data.addMulti<int>("label_Z_cc");
  data.addMulti<int>("label_Z_qq");

  data.addMulti<int>("label_H_bb");
  data.addMulti<int>("label_H_cc");
  data.addMulti<int>("label_H_qqqq");

  data.addMulti<int>("label_QCD_bb");
  data.addMulti<int>("label_QCD_cc");
  data.addMulti<int>("label_QCD_b");
  data.addMulti<int>("label_QCD_c");
  data.addMulti<int>("label_QCD_others");

  data.addMulti<int>("sample_isQCD");

  // legacy labels
  data.addMulti<int>("fj_labelLegacy");

  // JMAR label
  data.addMulti<int>("fj_labelJMAR");

  // gen-matched particle (top/W/etc.)
  data.addMulti<float>("fj_gen_pt");
  data.addMulti<float>("fj_gen_eta");

  // fatjet kinematics
  data.addMulti<float>("fj_pt");
  data.addMulti<float>("fj_eta");
  data.addMulti<float>("fj_phi");
  data.addMulti<float>("fj_mass");

  // substructure
  data.addMulti<float>("fj_tau1");
  data.addMulti<float>("fj_tau2");
  data.addMulti<float>("fj_tau3");
  data.addMulti<float>("fj_tau21");
  data.addMulti<float>("fj_tau32");

  // soft drop
  data.addMulti<float>("fj_sdmass");

  // subjets: soft drop gives up to 2 subjets
  data.addMulti<float>("fj_n_sdsubjets");

  data.addMulti<float>("fj_sdsj1_pt");
  data.addMulti<float>("fj_sdsj1_eta");
  data.addMulti<float>("fj_sdsj1_phi");
  data.addMulti<float>("fj_sdsj1_mass");
  data.addMulti<float>("fj_sdsj1_csv");
  data.addMulti<float>("fj_sdsj1_ptD");
  data.addMulti<float>("fj_sdsj1_axis1");
  data.addMulti<float>("fj_sdsj1_axis2");
  data.addMulti<float>("fj_sdsj1_mult");

  data.addMulti<float>("fj_sdsj2_pt");
  data.addMulti<float>("fj_sdsj2_eta");
  data.addMulti<float>("fj_sdsj2_phi");
  data.addMulti<float>("fj_sdsj2_mass");
  data.addMulti<float>("fj_sdsj2_csv");
  data.addMulti<float>("fj_sdsj2_ptD");
  data.addMulti<float>("fj_sdsj2_axis1");
  data.addMulti<float>("fj_sdsj2_axis2");
  data.addMulti<float>("fj_sdsj2_mult");

  // some variables used in a baseline tagger
  data.addMulti<float>("fj_ptDR");
  data.addMulti<float>("fj_relptdiff");
  data.addMulti<float>("fj_sdn2");


  //double-b
  data.addMulti<float>("fj_doubleb");

  //flavor info
  data.addMulti<int>("fj_isBB");
  data.addMulti<int>("fj_isNonBB");
  data.addMulti<int>("fj_nbHadrons");
  data.addMulti<int>("fj_ncHadrons");

  //double-b inputs
  data.addMulti<float>("fj_z_ratio");
  data.addMulti<float>("fj_trackSipdSig_3");
  data.addMulti<float>("fj_trackSipdSig_2");
  data.addMulti<float>("fj_trackSipdSig_1");
  data.addMulti<float>("fj_trackSipdSig_0");
  data.addMulti<float>("fj_trackSipdSig_1_0");
  data.addMulti<float>("fj_trackSipdSig_0_0");
  data.addMulti<float>("fj_trackSipdSig_1_1");
  data.addMulti<float>("fj_trackSipdSig_0_1");
  data.addMulti<float>("fj_trackSip2dSigAboveCharm_0");
  data.addMulti<float>("fj_trackSip2dSigAboveBottom_0");
  data.addMulti<float>("fj_trackSip2dSigAboveBottom_1");
  data.addMulti<float>("fj_tau1_trackEtaRel_0");
  data.addMulti<float>("fj_tau1_trackEtaRel_1");
  data.addMulti<float>("fj_tau1_trackEtaRel_2");
  data.addMulti<float>("fj_tau0_trackEtaRel_0");
  data.addMulti<float>("fj_tau0_trackEtaRel_1");
  data.addMulti<float>("fj_tau0_trackEtaRel_2");
  data.addMulti<float>("fj_tau_vertexMass_0");
  data.addMulti<float>("fj_tau_vertexEnergyRatio_0");
  data.addMulti<float>("fj_tau_vertexDeltaR_0");
  data.addMulti<float>("fj_tau_flightDistance2dSig_0");
  data.addMulti<float>("fj_tau_vertexMass_1");
  data.addMulti<float>("fj_tau_vertexEnergyRatio_1");
  data.addMulti<float>("fj_tau_flightDistance2dSig_1");
  data.addMulti<float>("fj_jetNTracks");
  data.addMulti<float>("fj_nSV");


}

bool FatJetInfoFiller::fill() {
  for (unsigned jetidx=0; jetidx<jets->size(); ++jetidx){
    const auto jet = jets->at(jetidx).correctedJet("Uncorrected"); // undo the JECs
    JetHelper jet_helper(&jet);
    // legacy label
    data.fillMulti<int>("fj_labelLegacy", fjmatch_.flavor(&jet, *genParticlesHandle).first);

    // JMAR label
    data.fillMulti<int>("fj_labelJMAR", fjmatch_.flavorJMAR(&jet, *genParticlesHandle, 0.6).first);


    // ----------------------------------------------------------------
    auto fjlabel = fjmatch_.flavorLabel(&jet, *genParticlesHandle, 0.6);

    data.fillMulti<int>("fj_label", fjlabel.first);

    data.fillMulti<int>("fj_isTop", fjlabel.first >= FatJetMatching::Top_all && fjlabel.first < FatJetMatching::W_all);
    data.fillMulti<int>("fj_isW",   fjlabel.first >= FatJetMatching::W_all && fjlabel.first < FatJetMatching::Z_all);
    data.fillMulti<int>("fj_isZ",   fjlabel.first >= FatJetMatching::Z_all && fjlabel.first < FatJetMatching::H_all);
    data.fillMulti<int>("fj_isH",   fjlabel.first >= FatJetMatching::H_all && fjlabel.first < FatJetMatching::QCD_all);
    data.fillMulti<int>("fj_isQCD", fjlabel.first >= FatJetMatching::QCD_all);

    data.fillMulti<int>("label_Top_bcq", fjlabel.first == FatJetMatching::Top_bcq);
    data.fillMulti<int>("label_Top_bqq", fjlabel.first == FatJetMatching::Top_bqq);
    data.fillMulti<int>("label_Top_bc",  fjlabel.first == FatJetMatching::Top_bc);
    data.fillMulti<int>("label_Top_bq",  fjlabel.first == FatJetMatching::Top_bq);

    data.fillMulti<int>("label_W_cq",    fjlabel.first == FatJetMatching::W_cq);
    data.fillMulti<int>("label_W_qq",    fjlabel.first == FatJetMatching::W_qq);

    data.fillMulti<int>("label_Z_bb",    fjlabel.first == FatJetMatching::Z_bb);
    data.fillMulti<int>("label_Z_cc",    fjlabel.first == FatJetMatching::Z_cc);
    data.fillMulti<int>("label_Z_qq",    fjlabel.first == FatJetMatching::Z_qq);

    data.fillMulti<int>("label_H_bb",    fjlabel.first == FatJetMatching::H_bb);
    data.fillMulti<int>("label_H_cc",    fjlabel.first == FatJetMatching::H_cc);
    data.fillMulti<int>("label_H_qqqq",  fjlabel.first == FatJetMatching::H_qqqq);

    data.fillMulti<int>("label_QCD_bb",  fjlabel.first == FatJetMatching::QCD_bb);
    data.fillMulti<int>("label_QCD_cc",  fjlabel.first == FatJetMatching::QCD_cc);
    data.fillMulti<int>("label_QCD_b",   fjlabel.first == FatJetMatching::QCD_b);
    data.fillMulti<int>("label_QCD_c",   fjlabel.first == FatJetMatching::QCD_c);
    data.fillMulti<int>("label_QCD_others", fjlabel.first == FatJetMatching::QCD_others);

    data.fillMulti<int>("sample_isQCD",  isQCDSample_);


    // gen-matched particle (top/W/etc.)
    data.fillMulti<float>("fj_gen_pt", fjlabel.second ? fjlabel.second->pt() : -999);
    data.fillMulti<float>("fj_gen_eta", fjlabel.second ? fjlabel.second->eta() : -999);
    // ----------------------------------------------------------------



    // fatjet kinematics
    data.fillMulti<float>("fj_pt", jet.pt());
    data.fillMulti<float>("fj_eta", jet.eta());
    data.fillMulti<float>("fj_phi", jet.phi());
    data.fillMulti<float>("fj_mass", jet.mass());

    // substructure
    //Possible userfloats (101X): NjettinessAK8Puppi:tau1 NjettinessAK8Puppi:tau2 NjettinessAK8Puppi:tau3 NjettinessAK8Puppi:tau4 ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1 ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2 ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3 ak8PFJetsCHSValueMap:NjettinessAK8CHSTau4 ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass ak8PFJetsCHSValueMap:eta ak8PFJetsCHSValueMap:mass ak8PFJetsCHSValueMap:phi ak8PFJetsCHSValueMap:pt ak8PFJetsPuppiSoftDropMass ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2 ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3 ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2 ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3 nb1AK8PuppiSoftDrop:ecfN2 nb1AK8PuppiSoftDrop:ecfN3 nb2AK8PuppiSoftDrop:ecfN2 nb2AK8PuppiSoftDrop:ecfN3 
    //Possible UserFloats (80X) are: NjettinessAK8:tau1 NjettinessAK8:tau2 NjettinessAK8:tau3 ak8PFJetsCHSPrunedMass ak8PFJetsCHSSoftDropMass ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1 ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2 ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3 ak8PFJetsPuppiValueMap:eta ak8PFJetsPuppiValueMap:mass ak8PFJetsPuppiValueMap:phi ak8PFJetsPuppiValueMap:pt 
    float tau1 = jet.userFloat("NjettinessAK8:tau1");
    float tau2 = jet.userFloat("NjettinessAK8:tau2");
    float tau3 = jet.userFloat("NjettinessAK8:tau3");
    data.fillMulti<float>("fj_tau1", tau1);
    data.fillMulti<float>("fj_tau2", tau2);
    data.fillMulti<float>("fj_tau3", tau3);
    data.fillMulti<float>("fj_tau21", tau1 > 0 ? tau2/tau1 : 1.01);
    data.fillMulti<float>("fj_tau32", tau2 > 0 ? tau3/tau2 : 1.01);

    // soft drop
    data.fillMulti<float>("fj_sdmass", jet.userFloat("ak8PFJetsCHSSoftDropMass"));

    // subjets: soft drop gives up to 2 subjets
    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("fj_n_sdsubjets", subjets.size());

    if (subjets.size() > 0){
      const auto &sj1 = subjets.at(0);
      JetHelper jh1(sj1);
      data.fillMulti<float>("fj_sdsj1_pt", sj1->pt());
      data.fillMulti<float>("fj_sdsj1_eta", sj1->eta());
      data.fillMulti<float>("fj_sdsj1_phi", sj1->phi());
      data.fillMulti<float>("fj_sdsj1_mass", sj1->mass());
      data.fillMulti<float>("fj_sdsj1_csv", sj1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      data.fillMulti<float>("fj_sdsj1_ptD", jh1.ptD());
      data.fillMulti<float>("fj_sdsj1_axis1", jh1.axis1());
      data.fillMulti<float>("fj_sdsj1_axis2", jh1.axis2());
      data.fillMulti<float>("fj_sdsj1_mult", jh1.mult());

      if (subjets.size() > 1){
        const auto &sj2 = subjets.at(1);
        JetHelper jh2(sj2);
        data.fillMulti<float>("fj_sdsj2_pt", sj2->pt());
        data.fillMulti<float>("fj_sdsj2_eta", sj2->eta());
        data.fillMulti<float>("fj_sdsj2_phi", sj2->phi());
        data.fillMulti<float>("fj_sdsj2_mass", sj2->mass());
        data.fillMulti<float>("fj_sdsj2_csv", sj2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        data.fillMulti<float>("fj_sdsj2_ptD", jh2.ptD());
        data.fillMulti<float>("fj_sdsj2_axis1", jh2.axis1());
        data.fillMulti<float>("fj_sdsj2_axis2", jh2.axis2());
        data.fillMulti<float>("fj_sdsj2_mult", jh2.mult());

        // some variables used in a baseline tagger
        float deltaR = reco::deltaR(*sj1, *sj2);
        float var_sd_0 = sj2->pt()/(sj1->pt()+sj2->pt());
        data.fillMulti<float>("fj_ptDR", jet.pt() * deltaR);
        data.fillMulti<float>("fj_relptdiff", std::abs(sj1->pt()-sj2->pt()) / jet.pt());
        data.fillMulti<float>("fj_sdn2", var_sd_0/std::pow(deltaR,-2));
      }
    }
  
    // --------
    // double-b

    const auto *bdsvTagInfo = jet.tagInfoBoostedDoubleSV(fjTagInfoName);
    assert(bdsvTagInfo);
    const auto &vars = bdsvTagInfo->taggingVariables();

    data.fillMulti<float>("fj_doubleb", jet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));

    //flavor info
    data.fillMulti<int>("fj_isBB", jet.jetFlavourInfo().getbHadrons().size() >= 2);
    data.fillMulti<int>("fj_isNonBB", jet.jetFlavourInfo().getbHadrons().size() < 2);
    data.fillMulti<int>("fj_nbHadrons", jet.jetFlavourInfo().getbHadrons().size());
    data.fillMulti<int>("fj_ncHadrons", jet.jetFlavourInfo().getcHadrons().size());

    //double-b inputs
    data.fillMulti<float>("fj_z_ratio", vars.get(reco::btau::z_ratio));
    data.fillMulti<float>("fj_trackSipdSig_3", vars.get(reco::btau::trackSip3dSig_3));
    data.fillMulti<float>("fj_trackSipdSig_2", vars.get(reco::btau::trackSip3dSig_2));
    data.fillMulti<float>("fj_trackSipdSig_1", vars.get(reco::btau::trackSip3dSig_1));
    data.fillMulti<float>("fj_trackSipdSig_0", vars.get(reco::btau::trackSip3dSig_0));
    data.fillMulti<float>("fj_trackSipdSig_1_0", vars.get(reco::btau::tau2_trackSip3dSig_0));
    data.fillMulti<float>("fj_trackSipdSig_0_0", vars.get(reco::btau::tau1_trackSip3dSig_0));
    data.fillMulti<float>("fj_trackSipdSig_1_1", vars.get(reco::btau::tau2_trackSip3dSig_1));
    data.fillMulti<float>("fj_trackSipdSig_0_1", vars.get(reco::btau::tau1_trackSip3dSig_1));
    data.fillMulti<float>("fj_trackSip2dSigAboveCharm_0", vars.get(reco::btau::trackSip2dSigAboveCharm));
    data.fillMulti<float>("fj_trackSip2dSigAboveBottom_0", vars.get(reco::btau::trackSip2dSigAboveBottom_0));
    data.fillMulti<float>("fj_trackSip2dSigAboveBottom_1", vars.get(reco::btau::trackSip2dSigAboveBottom_1));
    data.fillMulti<float>("fj_tau1_trackEtaRel_0", vars.get(reco::btau::tau2_trackEtaRel_0));
    data.fillMulti<float>("fj_tau1_trackEtaRel_1", vars.get(reco::btau::tau2_trackEtaRel_1));
    data.fillMulti<float>("fj_tau1_trackEtaRel_2", vars.get(reco::btau::tau2_trackEtaRel_2));
    data.fillMulti<float>("fj_tau0_trackEtaRel_0", vars.get(reco::btau::tau1_trackEtaRel_0));
    data.fillMulti<float>("fj_tau0_trackEtaRel_1", vars.get(reco::btau::tau1_trackEtaRel_1));
    data.fillMulti<float>("fj_tau0_trackEtaRel_2", vars.get(reco::btau::tau1_trackEtaRel_2));
    data.fillMulti<float>("fj_tau_vertexMass_0", vars.get(reco::btau::tau1_vertexMass));
    data.fillMulti<float>("fj_tau_vertexEnergyRatio_0", vars.get(reco::btau::tau1_vertexEnergyRatio));
    data.fillMulti<float>("fj_tau_vertexDeltaR_0", vars.get(reco::btau::tau1_vertexDeltaR));
    data.fillMulti<float>("fj_tau_flightDistance2dSig_0", vars.get(reco::btau::tau1_flightDistance2dSig));
    data.fillMulti<float>("fj_tau_vertexMass_1", vars.get(reco::btau::tau2_vertexMass));
    data.fillMulti<float>("fj_tau_vertexEnergyRatio_1", vars.get(reco::btau::tau2_vertexEnergyRatio));
    data.fillMulti<float>("fj_tau_flightDistance2dSig_1", vars.get(reco::btau::tau2_flightDistance2dSig));
    data.fillMulti<float>("fj_jetNTracks", vars.get(reco::btau::jetNTracks));
    data.fillMulti<float>("fj_nSV", vars.get(reco::btau::jetNSecondaryVertices));
  }
  return true;
}

} /* namespace deepntuples */


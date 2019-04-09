# HiggsToBBNtupleProducerTool

`ROOT` Ntuple producer for developing machine learning algorithms to differentiate Higgs to bb (Hbb) jets from quark or gluon jets (QCD).

## Setup

Running the data processing uses `CMSSW_8_0_28` 

```bash
cmsrel CMSSW_8_0_28
cd CMSSW_8_0_28/src/
cmsenv
git clone https://github.com/cms-legacydata-analyses/HiggsToBBNtupleProducerTool DeepNTuples -b opendata_80X
scram b -j 4
# if using CMS resources (like CRAB, you can set it up after:)
# voms-proxy-init -rfc -voms cms --valid 168:00
```

## Create `ROOT` ntuples

To run the python config, navigate to the `test` directory, and run it on an example file, by default `root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2016/RunIISummer16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0A83E4E2-34B6-E611-89A0-549F35AE4FA2.root`

```bash
cd DeepNTuples/NtupleAK8/test
cmsRun DeepNtuplizerAK8.py
```

Repeat this for all of the input files in the QCD and Hbb datasets.

## Merge outputs (with random mixing of different samples)

1. First create the file list (in this case we only take up to 3 output ROOT files for each QCD or Hbb sample).

```bash
cd DeepNTuples/NtupleAK8/run
export OUTDIR=/eos/uscms/store/group/lpcbtag/20181121_ak8_80x/OpenData
ls $OUTDIR > samples.txt
for i in `cat samples.txt`; 
 do cd ${OUTDIR}/${i}; 
 ls */*/*/output_*.root | head -3 > ${i}_max3files.txt;
 cd -; 
done
```

2. Merge the samples (with random mixing)

```bash
mergeSamples.py [events per output file] [output dir] [path to the filelist produced in step 1]
```
e.g.,
```bash
export MERGEDIR=/eos/uscms/store/group/lpcbtag/20181121_ak8_80x/merged_max3files
mergeSamples.py 200000 ${MERGEDIR} ${OUTDIR}/QCD_Pt_*/QCD_Pt_*max3files.txt ${OUTDIR}/Bulk*/Bulk*max3files.txt
```

3. Split into training and testing samples

```bash
export TRAINDIR=${MERGEDIR}/train
export TESTDIR=${MERGEDIR}/test
mkdir -p $TRAINDIR $TESTDIR
mv ${MERGEDIR}/ntuple_merged_[.0-8.].root ${TESTDIR}/
mv ${MERGEDIR}/ntuple_merged_*.root ${TRAINDIR}/
```

## Convert `ROOT` files to `HDF5` files using `uproot`

This step requires a more recent version of CMSSW.

```bash
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/
cmsenv
wget https://raw.githubusercontent.com/cms-legacydata-analyses/HiggsToBBNtupleProducerToo/opendata_80X/NtupleAK8/scripts/convert-uproot-opendata.py
```

Then you can run
```bash
python convert-uproot-opendata.py [input file (.root)] [output file (.h5)]
```
e.g.,
```
python convert-uproot-opendata.py ${TRAINDIR}/ntuple_merged_10.root ${TRAINDIR}/ntuple_merged_10.h5
```
which produces `HDF5` files with different arrays for each output variable. Note that during this conversion, only the information for up to 100 particle candidataes, 60 tracks, and 5 secondary vertices are saved in flattened, zero-padded, fixed-length arrays.

## Data set content

Variables are saved to the `ROOT` trees and `HDF5` tables jet by jet.

The reconstructed jets are clustered using the anti-kT algorithm with R=0.8 from particle flow (PF) candidates (AK8 jets). The standard L1+L2+L3+residual jet energy corrections are applied to the jets and pileup contamination is mitigated is using the charged hadron subtraction (CHS) algorithm. Selected features of inclusive (charged and neutral) PF candidates with pT > 0.95 GeV associated to the AK8 jet are provided. Additional selected features of charged PF candidates (i.e. tracks)  with pT > 0.95 GeV associated to the AK8 jet are also provided. Additional features of reconstructed secondary vertices associated to the AK8 jet (within deltaR < 0.8) are also provided.


| Data variable | Type | Description |
| :---------------------- | -----------------: | :---------------------- |
| event_no | UInt_t | Event number |
| npv | Float_t | Number of reconstructed primary vertices |
| ntrueInt | Float_t | True mean number of the poisson distribution for this event from which the number of interactions in each bunch crossing has been sampled |
| rho | Float_t | Median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates of the event |
| sample_isQCD | Float_t | Boolean that is 1 if the simulated sample corresponds to QCD multijet production |
| fj_doubleb | Float_t | Double-b tagging discriminant based on a boosted decision tree calculated for the AK8 jet (see [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/)) |
| fj_eta | Float_t | Pseudorapdity η of the AK8 jet |
| fj_gen_eta | Float_t | Pseudorapdity η of the generator-level, matched heavy particle: H, W, Z, top, etc. (default = -999)  |
| fj_gen_pt | Float_t | Transverse momentum of the generator-level, geometrically matched heavy particle: H, W, Z, t, etc. (default = -999)  |
| fj_isBB | Int_t | Boolean that is 1 if two or more b hadrons are clustered within the AK8 jet (see [SWGuideBTagMCTools](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)) |
| fj_isNonBB | Int_t | Boolean that is 1 if fewer than two b hadrons are clustered within the AK8 jet (see [SWGuideBTagMCTools](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)) |
| fj_nbHadrons | Int_t | Number of b hadrons that are clustered within the AK8 jet (see [SWGuideBTagMCTools](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)) |
| fj_ncHadrons | Int_t | Number of c hadrons that are clustered within the AK8 jet (see [SWGuideBTagMCTools](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)) |
| fj_isH | Int_t | Boolean that is 1 if a generator-level Higgs boson and its daughters are geometrically matched to the AK8 jet |
| fj_isTop | Int_t | Boolean that is 1 if a generator-level top quark and its daughters are geometrically matched to the AK8 jet |
| fj_isW | Int_t | Boolean that is 1 if a generator-level W boson and its daughters are geometrically matched to the AK8 jet |
| fj_isZ | Int_t | Boolean that is 1 if a generator-level Z boson and its daughters are geometrically matched to the AK8 jet |
| fj_isQCD | Int_t | Boolean that is 1 if none of the above matching criteria are satisified (H, top, W, Z) |
| fj_label | Int_t | Integer label: `Invalid=0,  Top_all=10, Top_bcq=11, Top_bqq=12, Top_bc=13, Top_bq=14, W_all=20, W_cq=21, W_qq=22, Z_all=30, Z_bb=31, Z_cc=32, Z_qq=33, H_all=40, H_bb=41, H_cc=42, H_qqqq=43, QCD_all=50, QCD_bb=51, QCD_cc=52, QCD_b=53, QCD_c=54, QCD_others=55` |
| fj_labelJMAR | Int_t | Alternative integer label from the CMS Jet/MET and Resolution (JMAR) group:  `Default=0, Top=1, W=2, Z=3, H=4` |
| fj_labelLegacy | Int_t | Alternative (legacy) integer label: `Default=0, Top=1, W=2, Z=3, H=4` |
| fj_jetNTracks | Float_t | Number of tracks associated with the AK8 jet |
| fj_nSV | Float_t | Number of secondary vertices associated with the AK8 jet (∆R < 0.7)|
| fj_n_sdsubjets | Float_t | Number of soft drop subjets in the AK8 jet (up to 2) |
| fj_mass | Float_t | Ungroomed mass of the AK8 jet |
| fj_phi | Float_t | Azimuthal angle ϕ of the AK8 jet |
| fj_pt | Float_t | Transverse momentum of the AK8 jet |
| fj_tau1 | Float_t | |
| fj_tau2 | Float_t | |
| fj_tau3 | Float_t | |
| fj_tau21 | Float_t | |
| fj_tau32 | Float_t | |
| fj_sdmass | Float_t | Soft drop mass of the AK8 jet |
| fj_ptDR | Float_t | Transverse momentum times the ΔR between the two soft drop subjets |
| fj_relptdiff | Float_t | Absolute relative difference between the transverse momenta of the two softdrop subjets |
| fj_sdn2 | Float_t | Fraction of second subjet transverse momentum times ∆R squared |
| fj_sdsj1_axis1 | Float_t | |
| fj_sdsj1_axis2 | Float_t | |
| fj_sdsj1_csv | Float_t | |
| fj_sdsj1_eta | Float_t | |
| fj_sdsj1_mass | Float_t | |
| fj_sdsj1_mult | Float_t | |
| fj_sdsj1_phi | Float_t | |
| fj_sdsj1_pt | Float_t | |
| fj_sdsj1_ptD | Float_t | |
| fj_sdsj2_axis1 | Float_t | |
| fj_sdsj2_axis2 | Float_t | |
| fj_sdsj2_csv | Float_t | |
| fj_sdsj2_eta | Float_t | |
| fj_sdsj2_mass | Float_t | |
| fj_sdsj2_mult | Float_t | |
| fj_sdsj2_phi | Float_t | |
| fj_sdsj2_pt | Float_t | |
| fj_sdsj2_ptD | Float_t | |
| fj_z_ratio | Float_t | z ratio varible as defined in [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/) |
| fj_trackSipdSig_0 | Float_t | First largest track 3D signed impact parameter significance (see [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/) ) |
| fj_trackSipdSig_1 | Float_t | Second largest track 3D signed impact parameter significance (see [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/) ) |
| fj_trackSipdSig_2 | Float_t | Third largest track 3D signed impact parameter significance (see [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/) )  |
| fj_trackSipdSig_3 | Float_t | Fourth largest track 3D signed impact parameter significance (see [CMS-BTV-16-002](http://cms-results.web.cern.ch/cms-results/public-results/publications/BTV-16-002/) )  |
| fj_trackSipdSig_0_0 | Float_t | First largest track 3D signed impact parameter significance associated to the first N-subjettiness axis |
| fj_trackSipdSig_0_1 | Float_t | Second largest track 3D signed impact parameter significance associated to the first N-subjettiness axis |
| fj_trackSipdSig_1_0 | Float_t | First largest track 3D signed impact parameter significance associated to the second N-subjettiness axis |
| fj_trackSipdSig_1_1 | Float_t | Second largest track 3D signed impact parameter significance associated to the second N-subjettiness axis |
| fj_trackSip2dSigAboveCharm_0 | Float_t | Track 2D signed impact parameter significance of first track lifting the combined invariant mass of the tracks above the c hadron threshold mass (1.5 GeV) |
| fj_trackSip2dSigAboveBottom_0 | Float_t | Track 2D signed impact parameter significance of first track lifting the combined invariant mass of the tracks above b hadron threshold mass (5.2 GeV) |
| fj_trackSip2dSigAboveBottom_1 | Float_t | Track 2D signed impact parameter significance of second track lifting the combined invariant mass of the tracks above b hadron threshold mass (5.2 GeV)  |
| fj_tau0_trackEtaRel_0 | Float_t | Smallest track pseudorapidity, relative to the jet axis, associated to the first N-subjettiness axis |
| fj_tau0_trackEtaRel_1 | Float_t | Second smallest track pseudorapidity, relative to the jet axis, associated to the first N-subjettiness axis  |
| fj_tau0_trackEtaRel_2 | Float_t | Third smallest track pseudorapidity, relative to the jet axis, associated to the first N-subjettiness axis  |
| fj_tau1_trackEtaRel_0 | Float_t | Smallest track pseudorapidity, relative to the jet axis, associated to the second N-subjettiness axis |
| fj_tau1_trackEtaRel_1 | Float_t | Second smallest track pseudorapidity, relative to the jet axis, associated to the second N-subjettiness axis |
| fj_tau1_trackEtaRel_2 | Float_t | Third smallest track pseudorapidity, relative to the jet axis, associated to the second N-subjettiness axis  |
| fj_tau_vertexMass_0 | Float_t | Total secondary vertex mass for the first N-subjettiness axis, defined as the invariant mass of all tracks from secondary vertices associated with the first N-subjettiness axis |
| fj_tau_vertexMass_1 | Float_t | Total secondary vertex mass for the second N-subjettiness axis, defined as the invariant mass of all tracks from secondary vertices associated with the second N-subjettiness axis |
| fj_tau_vertexEnergyRatio_0 | Float_t | Secondary vertex energy ratio for the first N-subjettiness axis, defined as the total energy of all secondary vertices associated with the first N-subjettiness axis divided by the total energy of all the tracks associated with the AK8 jet that are consistent with the primary vertex |
| fj_tau_vertexEnergyRatio_1 | Float_t | Secondary vertex energy ratio for the second N-subjettiness axis, defined as the total energy of all secondary vertices associated with the first N-subjettiness axis divided by the total energy of all the tracks associated with the AK8 jet that are consistent with the primary vertex |
| fj_tau_flightDistance2dSig_0 | Float_t | Transverse (2D) flight distance significance between the primary vertex and the secondary vertex with the smallest uncertainty on the 3D flight distance associated to the first N-subjettiness axis |
| fj_tau_flightDistance2dSig_1 | Float_t | Transverse (2D) flight distance significance between the primary vertex and the secondary vertex with the smallest uncertainty on the 3D flight distance associated to the second N-subjettiness axis |
| fj_tau_vertexDeltaR_0 | Float_t | Pseudoangular distance (∆R) between the first N-subjettiness axis and secondary vertex direction |
| n_pfcands | Int_t | Number of particle flow (PF) candidates associated to the AK8 jet with transverse momentum greater than 0.95 GeV |
| npfcands | Float_t | Number of particle flow (PF) candidates associated to the AK8 jet with transverse momentum greater than 0.95 GeV |
| pfcand_VTX_ass | Int_t | Primary vertex association quality for the PF candiate: `NotReconstructedPrimary=0, OtherDeltaZ=1, CompatibilityBTag=4, CompatibilityDz=5, UsedInFitLoose=6, UsedInFitTight=7`|
| pfcand_charge | Float_t | Electric charge of the PF candidate  |
| pfcand_deltaR | Float_t | Pseudoangular distance (∆R) between the PF candidate and the AK8 jet axis |
| pfcand_drminsv | Float_t | Minimum pseudoangular distance (∆R) between the secondary vertices and the PF candidate |
| pfcand_drsubjet1 | Float_t | Pseudoangular distance (∆R)  between the PF candidate and the first soft drop subjet |
| pfcand_drsubjet2 | Float_t | Pseudoangular distance (∆R)  between the PF candidate and the second soft drop subjet |
| pfcand_dxy | Float_t | Transverse impact paramater of the PF candidate, defined as the distance of closest approach of the PF candidate trajectory to the beam line in the transverse plane to the beam |
| pfcand_dxysig | Float_t | Transverse impact paramater significance of the PF candidate |
| pfcand_dz | Float_t | Longitudinal impact parameter, defined as the distance of closest approach of the PF candidate trajectory to the primary vertex projected on to the z direction |
| pfcand_dzsig | Float_t | Longitudinal impact parameter significance of the PF candidate |
| pfcand_erel | Float_t | Energy of the PF candidate divided by the energy of the AK8 jet |
| pfcand_etarel | Float_t | Pseudorapidity of the PF candidate relative to the AK8 jet axis |
| pfcand_phirel | Float_t | Pseudorapidity of the PF candidate relative to the pseudorapidity of the AK8 jet axis |
| pfcand_ptrel | Float_t | |
| pfcand_fromPV | Float_t | Integer indicating whether the PF candidate is consistent with the primary vertex: `NoPV=0, PVLoose=1, PVTight=2, PVUsedInFit=3` |
| pfcand_hcalFrac | Float_t | Fraction of energy of the PF candidate deposited in the hadron calorimeter |
| pfcand_isChargedHad | Float_t | Boolean that is 1 if the PF candidate is classified as a charged hadron |
| pfcand_isEl | Float_t | Boolean that is 1 if the PF candidate is classified as an electron |
| pfcand_isGamma | Float_t | Boolean that is 1 if the PF candidate is classified as an photon |
| pfcand_isMu | Float_t | Boolean that is 1 if the PF candidate is classified as an muon |
| pfcand_isNeutralHad | Float_t | Boolean that is 1 if the PF candidate is classified as a neutral hadron |
| pfcand_lostInnerHits | Float_t | Integer indicating whether any hits were lost in the inner silicon tracker: `validHitInFirstPixelBarrelLayer=-1, noLostInnerHits=0, oneLostInnerHit=1, moreLostInnerHits=2` |
| pfcand_mass | Float_t | Mass of the PF candidate |
| pfcand_puppiw | Float_t | Pileup per-particle identification (PUPPI) weight indicating whether the PF candidate is pileup-like (0) or not (1) |
| n_tracks | Int_t | Number of tracks associated with the AK8 jet |
| ntracks | Float_t | Number of tracks associated with the AK8 jet |
| trackBTag_DeltaR | Float_t | |
| trackBTag_Eta | Float_t | |
| trackBTag_EtaRel | Float_t | |
| trackBTag_JetDistVal | Float_t | |
| trackBTag_Momentum | Float_t | |
| trackBTag_PPar | Float_t | |
| trackBTag_PParRatio | Float_t | |
| trackBTag_PtRatio | Float_t | |
| trackBTag_PtRel | Float_t | |
| trackBTag_Sip2dSig | Float_t | |
| trackBTag_Sip2dVal | Float_t | |
| trackBTag_Sip3dSig | Float_t | |
| trackBTag_Sip3dVal | Float_t | |
| track_VTX_ass | Float_t | |
| track_charge | Float_t | |
| track_deltaR | Float_t | |
| track_detadeta | Float_t | |
| track_dlambdadz | Float_t | |
| track_dphidphi | Float_t | |
| track_dphidxy | Float_t | |
| track_dptdpt | Float_t | |
| track_drminsv | Float_t | |
| track_drsubjet1 | Float_t | |
| track_drsubjet2 | Float_t | |
| track_dxy | Float_t | |
| track_dxydxy | Float_t | |
| track_dxydz | Float_t | |
| track_dxysig | Float_t | |
| track_dz | Float_t | |
| track_dzdz | Float_t | |
| track_dzsig | Float_t | |
| track_erel | Float_t | |
| track_etarel | Float_t | |
| track_fromPV | Float_t | |
| track_isChargedHad | Float_t | |
| track_isEl | Float_t | |
| track_isMu | Float_t | |
| track_lostInnerHits | Float_t | |
| track_mass | Float_t | |
| track_normchi2 | Float_t | |
| track_phirel | Float_t | |
| track_pt | Float_t | |
| track_ptrel | Float_t | |
| track_puppiw | Float_t | |
| track_quality | Float_t | |
| n_sv | Int_t | Number of secondary vertices associated with the AK8 jet (∆R < 0.8) |
| nsv | Float_t | Number of secondary vertices associated with the AK8 jet (∆R < 0.8)|
| sv_chi2 | Float_t | |
| sv_costhetasvpv | Float_t | |
| sv_d3d | Float_t | |
| sv_d3derr | Float_t | |
| sv_d3dsig | Float_t | |
| sv_deltaR | Float_t | |
| sv_dxy | Float_t | |
| sv_dxyerr | Float_t | |
| sv_dxysig | Float_t | |
| sv_erel | Float_t | |
| sv_etarel | Float_t | |
| sv_mass | Float_t | |
| sv_ndf | Float_t | |
| sv_normchi2 | Float_t | |
| sv_ntracks | Float_t | |
| sv_phirel | Float_t | |
| sv_pt | Float_t | |
| sv_ptrel | Float_t | |
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


| Data variable | Type | Description |
| :---------------------- | -----------------: | :---------------------- |
| event_no | UInt_t | Event number |
| npv | Float_t | Number of reconstructed primary vertices |
| ntrueInt | Float_t | Number of true interactions |
| rho | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates of the event |
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
| fj_jetNTracks | Float_t | Number of tracks associated to the AK8 jet |
| fj_nSV | Float_t | Number of reconstructed possible secondary vertices in the AK8 jet |
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
| fj_sdn2 | Float_t | Fraction of second subjet transverse momentum times deltaR squared |
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
| fj_trackSipdSig_3 | Float_t | |
| fj_trackSipdSig_2 | Float_t | |
| fj_trackSipdSig_1 | Float_t | |
| fj_trackSipdSig_0 | Float_t | |
| fj_trackSipdSig_1_0 | Float_t | |
| fj_trackSipdSig_0_0 | Float_t | |
| fj_trackSipdSig_1_1 | Float_t | |
| fj_trackSipdSig_0_1 | Float_t | |
| fj_trackSip2dSigAboveCharm_0 | Float_t | |
| fj_trackSip2dSigAboveBottom_0 | Float_t | |
| fj_trackSip2dSigAboveBottom_1 | Float_t | |
| fj_tau1_trackEtaRel_0 | Float_t | |
| fj_tau1_trackEtaRel_1 | Float_t | |
| fj_tau1_trackEtaRel_2 | Float_t | |
| fj_tau0_trackEtaRel_0 | Float_t | |
| fj_tau0_trackEtaRel_1 | Float_t | |
| fj_tau0_trackEtaRel_2 | Float_t | |
| fj_tau_vertexMass_0 | Float_t | |
| fj_tau_vertexEnergyRatio_0 | Float_t | |
| fj_tau_vertexDeltaR_0 | Float_t | |
| fj_tau_flightDistance2dSig_0 | Float_t | |
| fj_tau_vertexMass_1 | Float_t | |
| fj_tau_vertexEnergyRatio_1 | Float_t | |
| fj_tau_flightDistance2dSig_1 | Float_t | |
| n_pfcands | Int_t | Number of particle flow candidates associated to the jet |
| npfcands | Float_t | Number of particle flow candidates associated to the jet |
| pfcand_VTX_ass | Float_t | |
| pfcand_charge | Float_t | |
| pfcand_deltaR | Float_t | |
| pfcand_drminsv | Float_t | |
| pfcand_drsubjet1 | Float_t | |
| pfcand_drsubjet2 | Float_t | |
| pfcand_dxy | Float_t | |
| pfcand_dxysig | Float_t | |
| pfcand_dz | Float_t | |
| pfcand_dzsig | Float_t | |
| pfcand_erel | Float_t | |
| pfcand_etarel | Float_t | |
| pfcand_fromPV | Float_t | |
| pfcand_hcalFrac | Float_t | |
| pfcand_isChargedHad | Float_t | |
| pfcand_isEl | Float_t | |
| pfcand_isGamma | Float_t | |
| pfcand_isMu | Float_t | |
| pfcand_isNeutralHad | Float_t | |
| pfcand_lostInnerHits | Float_t | |
| pfcand_mass | Float_t | |
| pfcand_phirel | Float_t | |
| pfcand_ptrel | Float_t | |
| pfcand_puppiw | Float_t | |
| n_tracks | Int_t | |
| ntracks | Float_t | |
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
| n_sv | Int_t | |
| nsv | Float_t | |
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
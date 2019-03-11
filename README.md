# DNNTuplesAK8

## Setup

Running the data processing uses `CMSSW_8_0_28` 

```
cmsrel CMSSW_8_0_28
cd CMSSW_8_0_28/src/
cmsenv
git clone https://github.com/cms-legacydata-analyses/DNNTuplesAK8 DeepNTuples -b opendata_80X
scram b 
# set up CRAB env; run it after cmsenv
voms-proxy-init -rfc -voms cms --valid 168:00
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
export TRAINDIR=/eos/uscms/store/group/lpcbtag/20181121_ak8_80x/merged_max3files/train
export TESTDIR=/eos/uscms/store/group/lpcbtag/20181121_ak8_80x/merged_max3files/test
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
wget https://raw.githubusercontent.com/cms-legacydata-analyses/DNNTuplesAK8/opendata_80X/NtupleAK8/scripts/convert-uproot-opendata.py
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

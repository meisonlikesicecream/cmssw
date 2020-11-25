# Setup instructions

Possibly only correct for this branch.

```
# Follow instructions from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_11_1_3

cmsrel CMSSW_11_1_3
cd CMSSW_11_1_3/src
cmsenv
git cms-init
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v3.1.9

git cms-addpkg L1Trigger/L1TCommon

# Merge branch with latest PhaseI jet+sums emulator
git cms-merge-topic -u EmyrClement:PhaseIJetSums_EmulatorImprovements_11_1_3

```

To rerun track trigger and L1T, the cmsDriver command on the twiki may not work.  This appears to work though:

```
cmsDriver.py step1 --conditions auto:phase2_realistic_T15 -n 2 --era Phase2C9 --eventcontent FEVTDEBUGHLT --runUnscheduled --mc -s RAW2DIGI,L1TrackTrigger,L1 --datatier FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,L1Trigger/Configuration/customisePhase2TTNoMC.customisePhase2TTNoMC,Configuration/DataProcessing/Utils.addMonitoring --geometry Extended2026D49 --fileout file:out.root --no_exec --nThreads 8 --python step1_L1_ProdLike.py  --filein "/store/mc/Phase2HLTTDRWinter20DIGI/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/50000/00850711-FB6B-E04A-A51A-C9C1F699E2F0.root"
```
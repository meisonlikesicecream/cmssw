# Setup instructions

Possibly only correct for this branch.

```
# Follow instructions from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#Recipe_AN1

cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_6_1_patch2
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.39.0
git cms-addpkg L1Trigger/L1TCommon

# Merge branch with latest PhaseI jet+sums emulator
git cms-merge-topic -u EmyrClement:PhaseIJetSums

```
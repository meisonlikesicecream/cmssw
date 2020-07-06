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

# Add Katie's branch
# This is Simone's dev branch with jets and sums, where all commits have been squashed
# There will likely be one conflict in corrector.cc - accept the changes in Katie's branch
git cms-merge-topic -u kwalkingshaw:sums_devel_squash

```

Emulator-firmware agreement for this branch:
MET: 94.8%, Jets: 94.5%, HT: 97.3%, MHT: 95.8%
This is based on the emulator with the config : test_EmulateTriggerChain_cfg.py
This runs the emulator that matches what has been implemented on a board (July 2020), and produces uncalibrated jets and sums from these uncalibrated jets.  The jet calibration has not yet been implemented in the firmware.

The config test_ComputeCalibrated7x7JetsAndSums_cfg.py runs the 7x7 jet finding algo, runs the calibration, and calculates the sums.  Does this for both PF and PUPPI candidates
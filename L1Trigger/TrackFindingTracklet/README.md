Very basic setup instructions:
```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src/
cmsenv

git cms-merge-topic EmyrClement:hybrid
svn co svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/trunkSimpleCode9 .

scramv1 b -j 8

cd L1Trigger/TrackFindingTracklet/test/
cmsRun L1TrackNtupleMaker_cfg.py
```

Here is a list of other manual changes you may have to make:
- Input files may or may not include the stub and cluster truth assocation, so you may need to modify the path in L1TrackNtupleMaker_cfg.py and add/remove the corresponding paths
- The choice of TMTT fitter is hardocded in interface/L1TGeomBase.hh
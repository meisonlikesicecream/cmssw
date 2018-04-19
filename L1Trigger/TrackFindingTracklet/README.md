Very basic setup instructions:
```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src/
cmsenv

git cms-merge-topic EmyrClement:hybrid_trackletFitter
svn co svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/branch/ejclemen_hybrid .
scramv1 b -j 8

cd TMTrackTrigger/TMTrackFinder/test/
cmsRun L1TrackNtupleMaker_cfg.py
```

Here is a list of other manual changes you may have to make:
- Input files may or may not include the stub and cluster truth assocation, so you may need to modify the path in L1TrackNtupleMaker_cfg.py and add/remove the corresponding paths

The branch of the tracklet code, EmyrClement:hybrid_trackletFitter, is based on skinnari:Tracklet_932.  The differences between the two can be seen here:
https://github.com/skinnari/cmssw/compare/skinnari:Tracklet_932...EmyrClement:hybrid_trackletFitter

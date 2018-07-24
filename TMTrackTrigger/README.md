# Setup instructions to just run

To checkout and run, or if your modifications won't need to be put into the central repository, do:

```
cmsrel CMSSW_9_3_8
cd CMSSW_9_3_8/src
cmsenv

git cms-init
git remote add -t portTMTT TMTT git@github.com:EmyrClement/cmssw.git
git fetch TMTT
git cms-merge-topic EmyrClement:portTMTT

scramv1 b -j 8

cd TMTrackTrigger/TMTrackFinder/test/
cmsRun tmtt_tf_analysis_cfg.py inputMC=../../MCsamples/937/RelVal/TTbar/PU200.txt Events=10
```
# Setup instructions for making modifications

If you plan to make changes, the procedure is very similar, except you do ```git cms-checkout-topic``` instead of ```git cms-merge-topic```:

```
cmsrel CMSSW_9_3_8
cd CMSSW_9_3_8/src
cmsenv

git cms-init
git remote add -t portTMTT TMTT git@github.com:EmyrClement/cmssw.git
git fetch TMTT
git cms-checkout-topic EmyrClement:portTMTT
scramv1 b -j 8

cd TMTrackTrigger/TMTrackFinder/test/
cmsRun tmtt_tf_analysis_cfg.py inputMC=../../MCsamples/937/RelVal/TTbar/PU200.txt Events=10
```

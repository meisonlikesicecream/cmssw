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
At this point, you should be on a local branch called portTMTT.  You can make modifications, and commit your changes to this branch.  

Below is a simple example of making modifications, pushing them to your remote repository, and making a pull request back to the repository you want your changes to end up in.  In this example, we will use the following as the "central" repository (something equivalent of trunkSimpleCode9 in svn) : https://github.com/EmyrClement/cmssw  Follow the link, and fork the repository to your own account.  You then need to add your newly forked repository as a remote repository in your local working area, which we will call ```origin```:
```
git remote add origin <url>
```
You can get the url by clicking on "Clone or download" on the webpage for YOUR repository, and will be something like ```git@github.com:<GitHubUsername>/cmssw.git```

Lets change branch to one called "myChanges":
```
git checkout -b myChanges
```
Modify some files:
```
echo "#Hello World" >> TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
```
Check the status of your modifications:
```
git status
```
Which should show something like:
```
# On branch myChanges
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
#
no changes added to commit (use "git add" and/or "git commit -a")
```
You can undo (revert) your changes (as explained in the message for ```git status```) with ```git checkout -- <file1> <file2> ...```
To see your modifications, you can do:
```
git diff TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
#or just
#git diff
```
Which should show something like:
```diff --git a/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py b/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
index 7b693b9..fc59e07 100644
--- a/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
+++ b/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
@@ -189,3 +189,4 @@ if options.outputDataset == 1:
   process.writeDataset.outputCommands.append('keep  *_TTAssociator*_TML1Tracks*_*')
 
   process.schedule = cms.Schedule(process.p, process.pa, process.pd)
+#Hello World
```
Add the files, and commit:
```
git add TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
git commit -m "Modification to tmtt_tf_analysis_cfg comfig file"
```
Doing ```git status``` will now show:
```
# On branch myChanges
nothing to commit, working directory clean
```
Now push these changes to your remote repository.
```
git push origin myChanges 
```
You are now ready to make a pull request back to the central repository.  Go the webpage for your remote repository, and you should see a box stating you have just pushed some changes to the myChanges branch, and gives the option to "Compare & pull request".  Make the pull request to merge your changes in the myChanges branch to EmyrClement/cmssw:portTMTT.  You shoulw end up with a pull request that looks like this : https://github.com/EmyrClement/cmssw/pull/1

'''
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src/
cmsenv

git cms-merge-topic EmyrClement:hybrid
svn co svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/trunkSimpleCode9 .

scramv1 b -j 8
'''
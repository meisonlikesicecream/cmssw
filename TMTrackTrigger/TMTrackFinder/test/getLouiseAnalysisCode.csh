#!/bin/tcsh
######################################################################################
# This gets the official (CMS L1 Track group approved) code from Louise Skinnari for 
# analysis L1 track objects (TTTrack objects).
# 
# This is an alternative to checking it out from git.
######################################################################################

# Name of directory where Louise's analysis code is.
#setenv louiseGitDir "https://raw.githubusercontent.com/skinnari/cmssw/Tracklet_932/L1Trigger/TrackFindingTracklet/test"

# Or use Emyr Clement's variant of it.
setenv louiseGitDir "https://raw.githubusercontent.com/EmyrClement/cmssw/Tracklet_93X_forTMTT/L1Trigger/TrackFindingTracklet/test"

curl -k ${louiseGitDir}/L1TrackNtupleMaker.cc     >! L1TrackNtupleMaker.cc
curl -k ${louiseGitDir}/L1TrackNtuplePlot.C       >! L1TrackNtuplePlot.C

# Add Buildfile.xml to test/ area to compile Louise Skinnari's code.

cat >! BuildFile.xml << 'EOF1'
<library   file="*.cc" name="L1trackLouiseCode">
  <flags   EDM_PLUGIN="1"/>
  <use name="TMTrackTrigger/TMTrackFinder"/>
  <!-- # Add no-misleading-indentation option to avoid warnings about bug in Boost library. -->
  <flags CXXFLAGS="-g -Wno-unused-variable -Wno-misleading-indentation"/>
</library>
'EOF1'

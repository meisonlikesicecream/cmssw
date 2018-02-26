############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")

GEOMETRY = "D17"

 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

if GEOMETRY == "D10": 
    print "using geometry " + GEOMETRY + " (flat)"
    process.load('Configuration.Geometry.GeometryExtended2023D10Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D10_cff')
elif GEOMETRY == "D17":
    print "using geometry " + GEOMETRY + " (tilted)"
    process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
elif GEOMETRY == "TkOnly":
    print "using geometry " + GEOMETRY + " (tilted)"
    process.load('L1Trigger.TrackTrigger.TkOnlyTiltedGeom_cff')
else:
    print "this is not a valid geometry!!!"

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

if GEOMETRY == "D10": 
    #D10 (flat barrel)
    Source_Files = cms.untracked.vstring(
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/044925F3-5F2E-E711-A92B-0CC47A7AB7A0.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/42543260-602E-E711-A41C-0025905A48C0.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/72D6B8B7-612E-E711-A727-0CC47A4D7692.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/786DD1A7-612E-E711-84EB-0025905B855E.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/9255C706-602E-E711-BDE0-0025905A609A.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/BE75B044-602E-E711-8DC3-0CC47A4D7692.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/EA054BA3-5F2E-E711-B7F4-0025905A6084.root",
)
elif GEOMETRY == "D17":
    #D17 (tilted barrel -- latest and greatest with T5 tracker, see: https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_0_pre2/Configuration/Geometry/README.md)
    Source_Files = cms.untracked.vstring(
    #"/store/user/ejclemen/L1TT/RelVal_932/WithTruthAssociation/TTbar/PU200/output_0.root",
    #"/store/user/ejclemen/L1TT/RelVal_932/WithTruthAssociation/TTbar/PU200/output_1.root",
    #"/store/user/ejclemen/L1TT/RelVal_932/WithTruthAssociation/TTbar/PU200/output_2.root",
    "/store/user/skinnari/TTbar_93X/PU200/output_0.root"
)
elif GEOMETRY == "TkOnly":
    Source_Files = cms.untracked.vstring(
    "file:MuMinus_1to50_TkOnly.root",
    "file:MuPlus_1to50_TkOnly.root"
)

process.source = cms.Source("PoolSource", fileNames = Source_Files, skipEvents = cms.untracked.uint32(0))

process.TFileService = cms.Service("TFileService", fileName = cms.string('TTbar_'+GEOMETRY+'_PU200.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# L1 tracking
############################################################

# remake stubs ?
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

if GEOMETRY == "D10": 
    TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(False)

if GEOMETRY != "TkOnly":
    from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
    TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
if GEOMETRY == "D10": 
    TTTracksFromTracklet.trackerGeometry = cms.untracked.string("flat")
#TTTracksFromTracklet.asciiFileName = cms.untracked.string("evlist.txt")

process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
process.TTTracks = cms.Path(process.L1TrackletTracks*process.TrackTriggerAssociatorTracks)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       ## tracking in jets stuff (--> requires AK4 genjet collection present!)
                                       TrackingInJets = cms.bool(True),
                                       GenJetInputTag = cms.InputTag("ak4GenJets", ""),
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

# use this if you want to re-run the stub making
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksWithTruth,process.ana)

# use this if cluster/stub associators not available 
#process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksWithTruth,process.ana)

process.schedule = cms.Schedule(process.TTTracks,process.ana)
# process.schedule = cms.Schedule(process.TTTracksWithTruth,process.ana)


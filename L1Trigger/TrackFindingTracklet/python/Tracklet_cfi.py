import FWCore.ParameterSet.Config as cms
from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer_params

TTTracksFromTracklet = cms.EDProducer("L1TrackProducer",
                                      SimTrackSource = cms.InputTag("g4SimHits"),
                                      SimVertexSource = cms.InputTag("g4SimHits"),
                                      TTStubSource = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                      MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                      MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                      TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                      TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                      BeamSpotSource = cms.InputTag("offlineBeamSpot"),
                                      asciiFileName = cms.untracked.string(""),
                                      trackerGeometryType  = cms.untracked.string(""),  #tilted barrel is assumed, use "flat" if running on flat
                                      tmttSettings = TMTrackProducer_params # Include TMTT settings
    )

# Remove digitisation in TMTT modules
TMTrackProducer_params.StubDigitize.EnableDigitize  = cms.bool(False)

# Inflate residual cuts in TMTT fitters (whilst getting the chain to work)
TMTrackProducer_params.TrackFitSettings.ResidualCut = cms.double(20.0)
TMTrackProducer_params.TrackFitSettings.GeneralResidualCut = cms.double(20.0)
TMTrackProducer_params.TrackFitSettings.KillingResidualCut = cms.double(20.0)

# Tracklet goes down to 2 GeV by default
TMTrackProducer_params.HTArraySpecRphi.HoughMinPt         = cms.double(2.0)


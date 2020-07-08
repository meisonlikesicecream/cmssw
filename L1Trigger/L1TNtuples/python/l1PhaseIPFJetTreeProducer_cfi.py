import FWCore.ParameterSet.Config as cms

l1PhaseIPFJetTree = cms.EDAnalyzer("L1PhaseIPFJetTreeProducer",
   genJetToken     = cms.untracked.InputTag("ak4GenJetsNoNu"),
   genMetToken     = cms.untracked.InputTag("genMetTrue"),
   ak4L1PF = cms.untracked.InputTag("ak4PFL1PuppiCorrected"),
   pfMetToken = cms.untracked.InputTag("l1PFMetPuppi"),
   l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates"),
   l1PhaseIPFSums = cms.untracked.InputTag("Phase1L1TSumsProducer", "Sums"),
#   ak4L1PF = cms.InputTag("L1TCorrectedPFJetProducer", "ak4PFL1PuppiCorrected"),
   maxL1Extra = cms.uint32(20)
)

runmenutree=cms.Path(l1PhaseIPFJetTree)





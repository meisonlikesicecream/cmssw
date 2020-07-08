import FWCore.ParameterSet.Config as cms
from math import pi

from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer
from L1Trigger.L1CaloTrigger.Phase1L1TJetCalibrator_9x9Jets_cfi import Phase1L1TJetCalibrator as Phase1L1TJetCalibrator9x9
from L1Trigger.L1CaloTrigger.Phase1L1TSumsProducer_cfi import Phase1L1TSumsProducer

Phase1L1TJetProducer9x9 = Phase1L1TJetProducer.clone(
	  jetIEtaSize = cms.uint32(9),
	  jetIPhiSize = cms.uint32(9),
	  outputCollectionName = cms.string("UncalibratedPhase1L1TJetFromPfCandidates")
)

Phase1L1TJetCalibrator9x9.inputCollectionTag = cms.InputTag("Phase1L1TJetProducer9x9", "UncalibratedPhase1L1TJetFromPfCandidates", "")
Phase1L1TJetCalibrator9x9.outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")

Phase1L1TSumsProducer9x9 = Phase1L1TSumsProducer.clone(
  jetCollectionTag = cms.InputTag("Phase1L1TJetCalibrator9x9", "Phase1L1TJetFromPfCandidates"),
)

Phase1L1TJetsSequence9x9 = cms.Sequence(
  Phase1L1TJetProducer9x9 +
  Phase1L1TJetCalibrator9x9 + 
  Phase1L1TSumsProducer9x9
)
import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils # ADDED
from L1Trigger.L1CaloTrigger.Phase1L1TSumsProducer_cfi import Phase1L1TSumsProducer


process = cms.Process("MET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service('TFileService', fileName = cms.string("Sums_Calibrated7x7Jets_Histograms_TTBar_PU200_104XMTD_HFCut.root"))

fileList = FileUtils.loadListFromFile('ZGammaToInvisiblePt120_PU200.list')
readFiles = cms.untracked.vstring(*fileList)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:TTBAR_Merged.root"),
  # skipEvents = cms.untracked.uint32(150000)
  # skipEvents = cms.untracked.uint32(158423)
  # fileNames = cms.untracked.vstring("file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/TTBar_200_10_0_4_MTD/inputs104X_1.root"),
  #fileNames = cms.untracked.vstring(
  #  "file:pf500.root",
  #)
  # fileNames = readFiles
)

process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')


process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('ComputeJetsAndSums_TTBar_PU200_104XMTD_HFCut.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer_*_*",
    "keep *_Phase1PFL1TJetProducer_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetCalibrator_*_*",
    "keep *_Phase1L1TSumsProducer_*_*",
    "keep *_Phase1PFL1TSumsProducer_*_*",
  ),
)

process.SaveGenL1TSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
  genMETCollectionTag = cms.InputTag("genMetTrue"),
  l1tMETCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"),
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"),
  l1tHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums")
)

process.SaveGenL1TPFSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
  genMETCollectionTag = cms.InputTag("genMetTrue"),
  l1tMETCollectionTag = cms.InputTag("Phase1PFL1TSumsProducer", "Sums"),
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"),
  l1tHTCollectionTag = cms.InputTag("Phase1PFL1TSumsProducer", "Sums")
)

process.SaveL1TSums = cms.EDAnalyzer("SaveL1Sums",
  l1tMETCollectionTag1 = cms.InputTag("Phase1L1TSumsProducer", "Sums"),
  l1tMETCollectionTag2 = cms.InputTag("Phase1PFL1TSumsProducer", "Sums"),
  l1tMETBranchName1 = cms.string("puppiL1TMET"),
  l1tMETBranchName2 = cms.string("pfL1TMET"),
)

# process.Phase1L1TSumsProducer = Phase1L1TSumsProducer

# process.p = cms.Path(process.Phase1L1TSumsProducer)

# process.e = cms.EndPath(process.out)


# process.out = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string('myOutputFile.root'),
#   outputCommands = cms.untracked.vstring(
#     "drop *",
#     "keep *_Phase1L1TJetProducer_*_*",
#     "keep *_ak4GenJetsNoNu_*_*",
#     "keep *_Phase1L1TJetCalibrator_*_*",
#   ),
# )

# from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer


# process.p = cms.Path(process.Phase1L1TJetsSequence)

process.Phase1PFL1TJetProducer = process.Phase1L1TJetProducer.clone(inputCollectionTag = cms.InputTag("l1pfCandidates", "PF"))
process.Phase1PFL1TSumsProducer = process.Phase1L1TSumsProducer.clone(particleCollectionTag = cms.InputTag("l1pfCandidates", "PF"))

process.p = cms.Path(process.Phase1L1TJetsSequence + process.Phase1PFL1TJetProducer + process.Phase1PFL1TSumsProducer + process.SaveGenL1TSums + process.SaveGenL1TPFSums + process.SaveL1TSums)

process.e = cms.EndPath(process.out)
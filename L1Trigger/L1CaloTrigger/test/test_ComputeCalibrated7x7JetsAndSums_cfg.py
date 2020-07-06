import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils # ADDED
from L1Trigger.L1CaloTrigger.Phase1L1TSumsProducer_cfi import Phase1L1TSumsProducer


process = cms.Process("MET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

# process.TFileService = cms.Service('TFileService', fileName = cms.string("Sums_Calibrated7x7Jets_Histograms_TTBar_PU200_104XMTD_HFCut.root"))

fileList = FileUtils.loadListFromFile('inputMCProcessed/TTbar_200PU_PhaseIITDRSpring19DR.txt')
readFiles = cms.untracked.vstring(*fileList)

process.source = cms.Source("PoolSource",
  # fileNames = cms.untracked.vstring("file:TTBAR_Merged.root"),
  # fileNames = cms.untracked.vstring('file:step2_2ev_reprocess_slim.root'),
  # fileNames = cms.untracked.vstring(
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_1.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_2.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_3.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_4.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_5.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_6.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_7.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_8.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_9.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_10.root',
  #   'file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_106X/TTbar_200PU/step2_2ev_reprocess_slim_11.root',
  #   ),
  fileNames = readFiles
  # skipEvents = cms.untracked.uint32(150000)
  # skipEvents = cms.untracked.uint32(158423)
  # fileNames = cms.untracked.vstring("file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/TTBar_200_10_0_4_MTD/inputs104X_1.root"),
  #fileNames = cms.untracked.vstring(
  #  "file:pf500.root",
  #)
  # fileNames = readFiles
   # inputCommands=cms.untracked.vstring(
   #     'keep *',
   #     'drop DetIdHGCSampleHGCDataFramesSorted_simHGCalUnsuppressedDigis_EE_*',
   #     # 'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusteredmNewDetSetVector_TTClustersFromPhase2TrackerDigis_ClusterInclusive_*',
   #     # 'drop l1tHGCalClusterBXVector_hgcalBackEndLayer1Producer_HGCalBackendLayer1Processor2DClustering_*',
   #     # 'drop l1tHGCalClusterBXVector_hgcalBackEndLayer1Producer_HGCalBackendLayer1Processor2DClustering_*',
   #     'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusterAssociationMap_TTClusterAssociatorFromPixelDigis_ClusterInclusive_*',
   #     # 'drop Phase2TrackerDigiedmDetSetVector_mix_Tracker_*',
   #     # 'drop l1tHGCalTriggerCellBXVector_hgcalConcentratorProducer_HGCalConcentratorProcessorSelection_*',
   #     # 'drop l1tHGCalTriggerCellBXVector_hgcalConcentratorProducer_HGCalConcentratorProcessorSelection_*',
   #     # 'drop EBDigiCollection_simEcalUnsuppressedDigis__*',
   #     # 'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusteredmNewDetSetVector_TTStubsFromPhase2TrackerDigis_ClusterAccepted_*',
   #     # 'drop DetIdHGCSampleHGCDataFramesSorted_simHGCalUnsuppressedDigis_HEfront_*',
   #     # 'drop FEDRawDataCollection_rawDataCollector__*',
   #     'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusterAssociationMap_TTClusterAssociatorFromPixelDigis_ClusterAccepted_*',
   #     # 'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubedmNewDetSetVector_TTStubsFromPhase2TrackerDigis_StubAccepted_*',
   #     # 'drop l1tHGCalMulticlusterBXVector_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_*',
   #     # 'drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis__*',
   #     # 'drop l1tHGCalMulticlusterBXVector_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_*',
   #     # 'drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis__*',
   #     # 'drop CSCDetIdCSCStripDigiMuonDigiCollection_simMuonCSCDigis_MuonCSCStripDigi_*',
   #     'drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubAssociationMap_TTStubAssociatorFromPixelDigis_StubAccepted_*',
   #     # 'drop GEMDigiSimLinkedmDetSetVector_simMuonME0Digis_ME0_*',
   #     # 'drop QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_*',
   #     # 'drop QIE11DataFrameHcalDataFrameContainer_simHcalDigis_HBHEQIE11DigiCollection_*',
   #     )
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

# process.SaveGenL1TSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("genMetTrue"),
#   l1tMETCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"),
#   genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"),
#   l1tHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums")
# )

# process.SaveGenL1TPFSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("genMetTrue"),
#   l1tMETCollectionTag = cms.InputTag("Phase1PFL1TSumsProducer", "Sums"),
#   genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"),
#   l1tHTCollectionTag = cms.InputTag("Phase1PFL1TSumsProducer", "Sums")
# )

# process.SaveL1TSums = cms.EDAnalyzer("SaveL1Sums",
#   l1tMETCollectionTag1 = cms.InputTag("Phase1L1TSumsProducer", "Sums"),
#   l1tMETCollectionTag2 = cms.InputTag("Phase1PFL1TSumsProducer", "Sums"),
#   l1tMETBranchName1 = cms.string("puppiL1TMET"),
#   l1tMETBranchName2 = cms.string("pfL1TMET"),
# )

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

# Add these to produce uncalibrated jets from PF (not PUPPI) inputs
# process.Phase1PFL1TJetProducer = process.Phase1L1TJetProducer.clone(inputCollectionTag = cms.InputTag("l1pfCandidates", "PF"))
# process.Phase1PFL1TSumsProducer = process.Phase1L1TSumsProducer.clone(particleCollectionTag = cms.InputTag("l1pfCandidates", "PF"))

process.load("L1Trigger.L1TNtuples.l1PhaseIITreeProducer_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple.root')
)

process.p = cms.Path(process.Phase1L1TJetsSequence )
# process.p = cms.Path(process.Phase1L1TJetsSequence )

process.e = cms.EndPath(process.out)
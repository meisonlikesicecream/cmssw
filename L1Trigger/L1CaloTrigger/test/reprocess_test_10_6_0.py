# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('REPR',eras.Phase2C8_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#process.Timing = cms.Service("Timing",
          #summaryOnly = cms.untracked.bool(False),
          #useJobReport = cms.untracked.bool(True)
#)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root'),
    #fileNames = cms.untracked.vstring('root://eoscms/eos/cms/store/relval/CMSSW_10_6_0_pre3/RelValMinBias_14TeV/GEN-SIM-DIGI-RAW/105X_upgrade2023_realistic_v5_2023D41noPU-v1/10000/E6CBA1C6-7A2E-A540-97B3-DE2C30AB70C8.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_pre3/RelValMuGunPt2To100/GEN-SIM-DIGI-RAW/105X_upgrade2023_realistic_v5_2023D41noPU-v2/10000/602E0B41-B698-6340-AC68-517578FEC457.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_106X_upgrade2023_realistic_v2_2023D41PU200-v1/10000/FEA5D564-937A-8D4B-9C9A-696EFC05AB58.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_6_0_patch2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/BC7B5A96-E3D2-ED48-81FC-35EF57134127.root'),
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/TTbar_14TeV_TuneCP5_Pythia8/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3_ext1-v3/60000/FFB5D0CA-208F-6040-A9BF-3F5354D0AA59.root'),
    # fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/Nu_E10-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v3/70000/0749C401-444F-E246-B1DF-134601C66363.root'),

    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step2_2ev_reprocess_slim.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.FEVTDEBUGHLToutput.outputCommands.append("drop *_*_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfCandidates_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_ak4PFL1PuppiCorrected_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TrackerEtMiss_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TrackerHtMiss_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_TwoLayerJets_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkCaloJets_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_Phase1L1TJetCalibrator_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1PFMetPuppi_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_ak4GenJetsNoNu_*_*")

process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1CaloJetProducer_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1CaloJetHTTProducer_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *PFTrack*_*_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *TTTrack*_*_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_*TTTrack*_*_*")
process.FEVTDEBUGHLToutput.outputCommands.append("keep *_addPileupInfo_*_*")

# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_simGmtStage2Digis_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1EGammaClusterEmuProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkElectronsCrystal_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkElectronsEllipticMatchCrystal_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkPhotonsCrystal_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1EGammaEEProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkElectronsHGC_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkElectronsEllipticMatchHGC_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkPhotonsHGC_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TrackerTaus_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkCaloTaus_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkEGTaus_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkGlbMuons_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkMuons_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TPSMuons_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1StubMatchedMuons_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkMuonStubEndCapS12_*_*")

# # process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubOverlap_*_*")
# # process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkMuonStubOverlap_MuonTracks_*")
# # process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkMuonStubOverlap_HscpTracks_*")

# process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubEndCap_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubOverlap_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_TwoLayerJets_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkCaloJets_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TrackerEtMiss_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TrackerHTMiss_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_ak4PFL1PuppiCorrected_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_ak4PFL1PuppiForMETCorrected_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_Phase1L1TJetProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_Phase1L1TJetCalibrator_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfCandidates_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfCandidates_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfSumsProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1CaloJetProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1CaloJetProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1CaloJetHTTProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_simKBmtfDigis_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_simOmtfDigis_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_simEmtfDigis_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1PFMetPuppi_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfProducerBarrel_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_VertexProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_VertexProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1TkPrimaryVertex_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1pfTauProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1NNTauProducerPuppi_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1NNTauProducer_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_L1HPSPFTauProducerPF_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkBsCandidates_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkBsCandidatesLooseWP_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("keep *_l1TkBsCandidatesTightWP_*_*")


# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_L1*_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_l1*_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubOverlap_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubEndCap_*_*")
# process.FEVTDEBUGHLToutput.outputCommands.append("drop *_l1TkMuonStubOverlap_*_*")
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TwoLayerJets_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_ak4PFL1PuppiCorrected_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_Phase1L1T_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_simK*mtfDigis_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('keep *_VertexProducer_*_*')
# # process.FEVTDEBUGHLToutput.outputCommands.append('drop *_*_*ayes*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_*_MergedTrackTruth_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiPixelDigis_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simHGCalUnsuppressedDigiss_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_hgcalVFEProducer_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_mix_MergedCaloTruth_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simHGCalUnsuppressedDigis_EE_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_hgcalBackEndLayer1Producer_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTClustersFromPhase2TrackerDigis_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTClusterAssociatorFromPixelDigis_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simEcalUnsuppressedDigis_*_*')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop *_hgcalConcentratorProducer_*_*')


# process.FEVTDEBUGHLToutput.outputCommands.append('drop Phase2TrackerDigiedmDetSetVector_mix_Tracker_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusteredmNewDetSetVector_TTStubsFromPhase2TrackerDigis_ClusterAccepted_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop DetIdHGCSampleHGCDataFramesSorted_simHGCalUnsuppressedDigis_HEfront_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop FEDRawDataCollection_rawDataCollector__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubedmNewDetSetVector_TTStubsFromPhase2TrackerDigis_StubAccepted_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tHGCalMulticlusterBXVector_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tHGCalMulticlusterBXVector_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis__REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop CSCDetIdCSCStripDigiMuonDigiCollection_simMuonCSCDigis_MuonCSCStripDigi_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubAssociationMap_TTStubAssociatorFromPixelDigis_StubAccepted_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop GEMDigiSimLinkedmDetSetVector_simMuonME0Digis_ME0_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop QIE11DataFrameHcalDataFrameContainer_simHcalDigis_HBHEQIE11DigiCollection_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop EBDigiCollection_simEcalDigis_ebDigis_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tHGCalTowerBXVector_hgcalTowerProducer_HGCalTowerProcessor_REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tHGCalTowerBXVector_hgcalTowerProducer_HGCalTowerProcessor_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0DigiPreRecoMuonDigiCollection_simMuonME0PseudoDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0RecHitsOwnedRangeMap_me0RecHits__REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop StripDigiSimLinkedmDetSetVector_simMuonME0Digis_ME0_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0RecHitsOwnedRangeMap_me0RecHitsCoarse__REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0RecHitsOwnedRangeMap_me0RecHitsCoarse__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop QIE10DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HFQIE10DigiCollection_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop QIE10DataFrameHcalDataFrameContainer_simHcalDigis_HFQIE10DigiCollection_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop HcalTriggerPrimitiveDigisSorted_simHcalTriggerPrimitiveDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop PileupSummaryInfos_addPileupInfo__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0DigiPreRecoMuonDigiCollection_simMuonME0PseudoReDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop HODataFramesSorted_simHcalUnsuppressedDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop EcalTriggerPrimitiveDigisSorted_simEcalTriggerPrimitiveDigis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0DigiPreRecoMuonDigiCollection_simMuonME0PseudoReDigisCoarse__REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop GlobalObjectMapRecord_hltGtStage2ObjectMap__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop GlobalObjectMapRecord_simGtStage2Digis__REPR')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop GlobalObjectMapRecord_simGtStage2Digis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0DigiPreRecoMuonDigiCollection_simMuonME0PseudoReDigisCoarse__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0DigiMuonDigiCollection_simMuonME0Digis__HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_MuonCSCStripDigiSimLinks_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_RPCDigiSimLink_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop GEMDigiSimLinkedmDetSetVector_simMuonGEMDigis_GEM_HLT')
# process.FEVTDEBUGHLToutput.outputCommands.append('drop ME0DetIdME0PadDigiMuonDigiCollection_simMuonME0PadDigis__HLT')

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.load("L1Trigger.L1TNtuples.l1PhaseIITreeProducer_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1NtuplePhaseII_160.root')
)



# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step,process.runmenutree,process.endjob_step,process.FEVTDEBUGHLToutput_step)
# process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet,L1TTurnOffHGCalTPs_v9,configureCSCLCTAsRun2
process = L1TrackTriggerTracklet(process)
#process = L1TTurnOffHGCalTPs_v9(process)
from L1Trigger.L1TMuonEndCap.customise_Phase2 import customise as customise_Phase2
process = customise_Phase2(process)


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

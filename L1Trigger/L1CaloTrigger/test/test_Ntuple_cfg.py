import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("Ntuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments

options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('outFile',
                 'L1Ntuple.root',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file')

options.parseArguments()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# fileList = FileUtils.loadListFromFile('snu.list')
# readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring( 'file:myOutputFile.root' )
)



process.load("L1Trigger.L1TNtuples.l1PhaseIPFJetTreeProducer_cfi")

process.l1PhaseIPFJetTree9x9 = process.l1PhaseIPFJetTree.clone(
   l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator9x9", "Phase1L1TJetFromPfCandidates"),
   l1PhaseIPFSums = cms.untracked.InputTag("Phase1L1TSumsProducer9x9", "Sums"),
   )


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple.root')
)

process.p = cms.Path(process.l1PhaseIPFJetTree+process.l1PhaseIPFJetTree9x9)

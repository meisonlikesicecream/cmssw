import FWCore.ParameterSet.Config as cms
from math import pi

caloEtaSegmentation = cms.vdouble(
  -5.0, -4.917, -4.833, -4.75, -4.667, -4.583, -4.5, -4.417, -4.333, -4.25, 
  -4.167, -4.083, -4.0, -3.917, -3.833, -3.75, -3.667, -3.583, -3.5, -3.417, 
  -3.333, -3.25, -3.167, -3.083, -3.0, -2.917, -2.833, -2.75, -2.667, -2.583, 
  -2.5, -2.417, -2.333, -2.25, -2.167, -2.083, -2.0, -1.917, -1.833, -1.75, 
  -1.667, -1.583, -1.5, -1.417, -1.333, -1.25, -1.167, -1.083, -1.0, -0.917, 
  -0.833, -0.75, -0.667, -0.583, -0.5, -0.417, -0.333, -0.25, -0.167, -0.083, 
  0.0, 0.083, 0.167, 0.25, 0.333, 0.417, 0.5, 0.583, 0.667, 0.75, 0.833, 0.917, 
  1.0, 1.083, 1.167, 1.25, 1.333, 1.417, 1.5, 1.583, 1.667, 1.75, 1.833, 1.917, 
  2.0, 2.083, 2.167, 2.25, 2.333, 2.417, 2.5, 2.583, 2.667, 2.75, 2.833, 2.917, 
  3.0, 3.083, 3.167, 3.25, 3.333, 3.417, 3.5, 3.583, 3.667, 3.75, 3.833, 3.917, 
  4.0, 4.083, 4.167, 4.25, 4.333, 4.417, 4.5, 4.583, 4.667, 4.75, 4.833, 4.917, 5.0)

#caloEtaSegmentation = cms.vdouble(
#    0.0, 0.0833, 0.1666, 0.2499, 0.3332, 0.4165, 0.4998, 0.5831, 0.6664, 0.7497, 
#    0.833, 0.9163, 0.9996, 1.0829, 1.1662, 1.2495, 1.3328, 1.4161, 1.5)

# lut configuration, you can generate your own with test/generateConSinPhi.py 
sinPhi = cms.vdouble(0.04374, 0.13087, 0.21701, 0.30149, 0.38365, 0.46289, 0.53858, 0.61015)
cosPhi = cms.vdouble(0.99904, 0.9914, 0.97617, 0.95347, 0.92348, 0.88642, 0.84257, 0.79229)

# sinPhi = cms.vdouble(0.046875, 0.132813, 0.21875, 0.304688, 0.382813, 0.460938, 0.539063, 0.609375 )
# cosPhi = cms.vdouble( 1, 0.992188, 0.976563, 0.953125, 0.921875, 0.882813, 0.84375, 0.789063 )


Phase1L1TJetProducer = cms.EDProducer('Phase1L1TJetProducer',
  inputCollectionTag = cms.InputTag("l1pfCandidates", "Puppi"),
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.15),
  phiUp = cms.double(3.15),
  # nBinsPhi = cms.uint32(8),
  # phiLow = cms.double(0),
  # phiUp = cms.double(0.7),
  jetIEtaSize = cms.uint32(7),
  jetIPhiSize = cms.uint32(7),
  trimmedGrid = cms.bool(False),
  seedPtThreshold = cms.double(5), # GeV
  puSubtraction = cms.bool(False),
  philsb = cms.double(0.0043633231),
  etalsb = cms.double(0.0043633231),
  outputCollectionName = cms.string("UncalibratedPhase1L1TJetFromPfCandidates"),
  vetoZeroPt = cms.bool(True),
  pfEtaRegions = cms.vdouble( -5., -4.5, -4., -3.5, -3., -2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5, 3., 3.5, 4., 4.5, 5. ),
  pfPhiRegions = cms.vdouble( -3.5, -2.8, -2.1, -1.4, -0.7, 0, 0.7, 1.4, 2.1, 2.8, 3.5 ),#, 4.2, 4.9, 5.6, 6.3 ),
  maxInputsPerPFRegion = cms.uint32( 24 ),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  metAbsEtaCut = cms.double(3),
  metHFAbsEtaCut = cms.double(5),
  debug = cms.bool(False)
)

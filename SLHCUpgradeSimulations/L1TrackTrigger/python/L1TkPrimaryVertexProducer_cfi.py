import FWCore.ParameterSet.Config as cms

L1TkPrimaryVertex = cms.EDProducer('L1TkFastVertexProducer',
     L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
     L1Tk_nPar = cms.int32( 4 ) ,   # use 4/5 parameter tracks?
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 1e10 ),
     PTMINTRA = cms.double( 2.),        # PTMIN of L1Tracks, in GeV
     nStubsmin = cms.int32( 4 ) ,       # minimum number of stubs
     nStubsPSmin = cms.int32( 3 ),       # minimum number of stubs in PS modules 
     nBinning = cms.int32( 601 ),        # number of bins for the temp histo (from -30 cm to + 30 cm)
     PTMAX = cms.double( 50. ),          # in GeV. When PTMAX > 0, tracks with PT above PTMAX are considered as
					 # mismeasured and are treated according to HighPtTracks below.
					 # When PTMAX < 0, no special treatment is done for high PT tracks.
					 # If PTMAX < 0, no saturation or truncation is done.
     HighPtTracks = cms.int32( 0 ),	 # when = 0 : truncation. Tracks with PT above PTMAX are ignored 
					 # when = 1 : saturation. Tracks with PT above PTMAX are set to PT=PTMAX.
     MonteCarloVertex = cms.bool( False ),    #  when True: dont run the vxt finding algo but pick up the MC generated vtx
     doPtComp = cms.bool( True ),       # track-stubs PT compatibility cut
     doTightChi2 = cms.bool( False )    # chi2dof < 5 for tracks with PT > 10
)
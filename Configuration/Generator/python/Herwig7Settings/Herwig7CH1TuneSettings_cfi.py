import FWCore.ParameterSet.Config as cms

herwig7CH1SettingsBlock = cms.PSet(
    herwig7CH1PDF = cms.vstring(
            'cd /Herwig/Partons',
            'create ThePEG::LHAPDF PDFSet_nnlo ThePEGLHAPDF.so',
            'set PDFSet_nnlo:PDFName NNPDF31_nnlo_as_0118.LHgrid',
            'set PDFSet_nnlo:RemnantHandler HadronRemnants',
            'set /Herwig/Particles/p+:PDF PDFSet_nnlo',
            'set /Herwig/Particles/pbar-:PDF PDFSet_nnlo',

            'set /Herwig/Partons/PPExtractor:FirstPDF  PDFSet_nnlo',
            'set /Herwig/Partons/PPExtractor:SecondPDF PDFSet_nnlo',

            'set /Herwig/Shower/ShowerHandler:PDFA PDFSet_nnlo',
            'set /Herwig/Shower/ShowerHandler:PDFB PDFSet_nnlo',
            
            'set /Herwig/Shower/ShowerHandler:PDFARemnant PDFSet_nnlo',
            'set /Herwig/Shower/ShowerHandler:PDFBRemnant PDFSet_nnlo',
            'set /Herwig/Partons/MPIExtractor:FirstPDF PDFSet_nnlo',
            'set /Herwig/Partons/MPIExtractor:SecondPDF PDFSet_nnlo', 

            'cd /',
        ),
    herwig7CH1AlphaS = cms.vstring(
        'cd /Herwig/Shower',
        'set AlphaQCD:AlphaMZ 0.118',
        'cd /'
        ),
    herwig7CH1MPISettings = cms.vstring(
        'read snippets/SoftModel.in',
        'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.4002',
        'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 2.322',
        'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.532',
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.1565',
        'set /Herwig/Partons/RemnantDecayer:ladderPower -0.08',
        'set /Herwig/Partons/RemnantDecayer:ladderNorm 0.95',
                                )
)

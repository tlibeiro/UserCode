import FWCore.ParameterSet.Config as cms

####Init Scales for clustering tree
fftjet_patreco_scales_25 = cms.PSet(
    Class = cms.string("EquidistantInLogSpace"),
    minScale = cms.double(0.087),
    maxScale = cms.double(0.224015352072),
    nScales = cms.uint32(25)
)
fftjet_patreco_scales_50 = cms.PSet(
    Class = cms.string("EquidistantInLogSpace"),
    minScale = cms.double(0.087),
    maxScale = cms.double(0.6),
    nScales = cms.uint32(50)
)

vplusganalyzer = cms.EDAnalyzer('VplusGAnalyzer'
)

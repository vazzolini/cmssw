import FWCore.ParameterSet.Config as cms

vertexCutsBlock = cms.PSet(
	vertexCuts = cms.PSet(
		fracPV = cms.double(0.65),
		distSig3dMax = cms.double(99999.9),
		distVal2dMax = cms.double(2.5),
		useTrackWeights = cms.bool(True),
		maxDeltaRToJetAxis = cms.double(0.5),
		v0Filter = cms.PSet(k0sMassWindow = cms.double(0.05)),
		distSig2dMin = cms.double(3.0),
		multiplicityMin = cms.uint32(2),
		massMax = cms.double(6.5),
		distSig2dMax = cms.double(99999.9),
		distVal3dMax = cms.double(99999.9),
		minimumTrackWeight = cms.double(0.5),
		distVal3dMin = cms.double(-99999.9),
		distVal2dMin = cms.double(0.01),
		distSig3dMin = cms.double(-99999.9)
	)
)

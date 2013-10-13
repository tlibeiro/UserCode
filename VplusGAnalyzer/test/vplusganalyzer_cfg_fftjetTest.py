import FWCore.ParameterSet.Config as cms
process = cms.Process("Analyzer")

wantSummary = True
numEventstoRun = -1
runOnMC = True
#outputFile = 'VplusGTree_fftjetTest_100Karis.root'
outputFile = 'OUTPUTFILE'

##----- Global tag: conditions database ------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
############################################
if not runOnMC:
    process.GlobalTag.globaltag = 'GR_R_53_V10::All'
else:
    process.GlobalTag.globaltag = 'START53_V7E::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numEventstoRun))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#									'file:/uscms/home/tlibeiro/VPlusGamma/CMSSW_5_3_8_patch1/src/TopQuarkAnalysis/TopPairBSM/test/tlbsm_53x_v3_mc_lvjj.root'
#'/store/user/smpjs/wgamma/ntran/GJets_HT-400ToInf_8TeV-madgraph/GJets_HT0400ToInf_w_53X_TLBSM_test2/2cb433dc57c751df38849560e5ef7270/tlbsm_53x_v3_mc_lvjj_63_1_mJX.root'
#/GJets_HT-400ToInf_8TeV-madgraph/StoreResults-Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3-99bd99199697666ff01397dad5652e9e/USER
#				'/store/user/lpctlbsm/meloam/GJets_HT-200To400_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_25_3_JpF.root',
#    '/store/user/smpjs/tlibeiro/GJets_HT-400ToInf_8TeV-madgraph/GJets_HT400ToInf_w_53X_TLBSM/82be9d6a1eead28db790c645493bb830/tlbsm_53x_v3_mc_lvjj_100_1_PHq.root'
###privately produced wa, waj, wajj samepl
#'/store/user/tlibeiro/WA_merged_w_GENSIM/WA_merged_w_53X_TLBSM/41d252faad44973494be0c526d69a069/tlbsm_53x_v3_mc_lvjj_11_1_mA1.root',
###signal from aris 100K 
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_1.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_2.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_3.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_4.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_5.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_6.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_7.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_8.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_9.root',
#'/store/user/tlibeiro/signalFromAris_100K/tlbsm_53x_v3_mc_lvjj_10.root'
##
'INPUTFILE'
		)
)
##---------  Load standard Reco modules ------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
##output file 
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outputFile),
    closeFileFast = cms.untracked.bool(True)
)
#############fftjet config####################################
from RecoJets.FFTJetProducers.fftjetcommon_cfi import *
calib_cone_radii = (0.2,0.3,0.4,0.5)
la_jets = (2,3,4)
#make new modules for each radius value
def adjust_cone_radius(R,pfJetsJetrecoPrototype,LAJets):
    #
    rtag = "R%d" % (R*100.0)
    latag = "LA"+str(LAJets)
    #
    newPFReco = pfJetsJetrecoPrototype.clone(
        recoScaleCalcPeak = cms.PSet(
            Class = cms.string("ConstDouble"),
            value = cms.double(R)
        ),
        recoScaleCalcJet = cms.PSet(
            Class = cms.string("ConstDouble"),
            value = cms.double(R)
        )
    )
    newPFReco_module = "fftPFJetReco"+latag+rtag
    #
    newPFReco.resolution = cms.string('locallyAdaptive')
    newPFReco.nClustersRequested = cms.uint32(LAJets)
    newPFCorr = configure_fftjet_correction_producer(
    (jetCorrectionSequenceTag,), newPFReco_module)
    newPFCorr.verbose = cms.untracked.bool(False)
    newPFCorr.calculatePileup = cms.bool(False) 
    newPFCorr.subtractPileup = cms.bool(False) 
    newPFCorr_module = "fftPFJetCorrected"+latag+rtag
 
    return ((newPFReco, newPFReco_module),
            (newPFCorr, newPFCorr_module)
            )

from TopQuarkAnalysis.MakeCandidateCollection.makecandidatecollection_cfi import * 
process.makeCandidates = makeCandidateCollection
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring( 'MakeCandidatesError'))
#############fftjet config####################################
# apply corrections
database = 'sqlite_file:fftjet_corr.db'
tableName = "L2L3"
tableCategory = "PF"
jetCorrectionSequenceTag = "PF0"
# Configure database access
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = database
# Configure event setup to work with jet corrections
from JetMETCorrections.FFTJetModules.fftjetcorrectionesproducer_cfi import *
configure_fftjet_pooldbessource(process, jetCorrectionSequenceTag)
config, esProducer = configure_L2L3_fftjet_esproducer(jetCorrectionSequenceTag,
                                                      tableName, tableCategory)
setattr(process, esProducer, config)
# Configure the jet corrections module
from RecoJets.FFTJetProducers.fftjetcorrectionproducer_cfi import *
#Using the table by FFTJetProducer:
#----------------------------------
from JetMETCorrections.FFTJetModules.fftjetlookuptableesproducer_cfi import *
##pucorrections 
etaSequenceTag = "EtaFlatteningFactors"
rhoSequenceTag = "PileupRhoCalibration"
configure_fftjetlut_pooldbessource(process, etaSequenceTag)
config, esProducer = configure_fftjetlut_esproducer(etaSequenceTag)
setattr(process, esProducer, config)
configure_fftjetlut_pooldbessource(process, rhoSequenceTag)
config, esProducer = configure_fftjetlut_esproducer(rhoSequenceTag)
setattr(process, esProducer, config)

fftjet_pileup_grid_pf_calibrated = cms.PSet(
    nEtaBins = cms.uint32(256),
    etaMin = cms.double(-math.pi*2),
    etaMax = cms.double(math.pi*2),
    nPhiBins = cms.uint32(128),
    phiBin0Edge = cms.double(0.0),
    title = cms.untracked.string("FFTJet Pileup Grid")
)
from RecoJets.FFTJetProducers.fftjetpileupestimator_pf_cfi import *
#Using the table by FFTJetPileupEstimator (on top of normal config):
#-------------------------------------------------------------------
fftjetPileupEstimatorPf.inputLabel = cms.InputTag("pileupprocessor", "FFTJetPileupPF")
fftjetPileupEstimatorPf.cdfvalue = cms.double(0.4)
fftjetPileupEstimatorPf.uncertaintyZones = cms.vdouble()
fftjetPileupEstimatorPf.calibrationCurve = cms.PSet(
    Class = cms.string("Polynomial"),
    c0 = cms.double(0.0)
)
fftjetPileupEstimatorPf.uncertaintyCurve = cms.PSet(
    Class = cms.string("Polynomial"),
    c0 = cms.double(0.0)
)
fftjetPileupEstimatorPf.calibTableRecord = cms.string(fftjet_lut_types[rhoSequenceTag].LUTRecord)
fftjetPileupEstimatorPf.calibTableCategory = cms.string("Pileup")
fftjetPileupEstimatorPf.uncertaintyZonesName = cms.string("FFTPileupRhoUncertaintyZones")
fftjetPileupEstimatorPf.calibrationCurveName = cms.string("FFTPileupRhoCalibrationTable")
fftjetPileupEstimatorPf.uncertaintyCurveName = cms.string("FFTPileupRhoUncertaintyTable")
fftjetPileupEstimatorPf.loadCalibFromDB = cms.bool(True)

from RecoJets.FFTJetProducers.fftjetpatrecoproducer_cfi import *
# Pattern recognition for PFJets
fftjetPFPatReco = fftjetPatrecoProducer.clone()
fftjetPFPatReco.src = cms.InputTag('makeCandidates:MyCandidates')
#fftjetPFPatReco.src = cms.InputTag('selectedPatJetsCA8PF:pfCandidates')
fftjetPFPatReco.jetType = cms.string('PFJet')
#fftjetPFPatReco.insertCompleteEvent = cms.bool(False)
fftjetPFPatReco.InitialScales = fftjet_patreco_scales_50
fftjetPFPatReco.calculateClusterRadii = cms.bool(True)
fftjetPFPatReco.calculateClusterSeparations = cms.bool(True)
# The Jet producer module
from RecoJets.FFTJetProducers.fftjetproducer_cfi import *
fftjetPFJetMaker = fftjetJetMaker.clone()
fftjetPFJetMaker.InitialScales = fftjet_patreco_scales_50
fftjetPFJetMaker.src = cms.InputTag('makeCandidates:MyCandidates')
fftjetPFJetMaker.jetType = cms.string('PFJet')
fftjetPFJetMaker.resolution = cms.string('fixed')
fftjetPFJetMaker.maxIterations = cms.uint32(1000)
fftjetPFJetMaker.fixedScale = cms.double(0.1)
fftjetPFJetMaker.PeakSelectorConfiguration = fftjet_peak_selector_allpass
fftjetPFJetMaker.treeLabel = cms.InputTag("fftjetPFPatReco", "FFTJetPatternRecognition")
fftjetPFJetMaker.calculatePileup = cms.bool(True)
fftjetPFJetMaker.subtractPileup = cms.bool(True)
##subtract pile up 
fftjetPFJetMaker.PileupGridConfiguration = fftjet_pileup_grid_pf_calibrated
fftjetPFJetMaker.pileupDensityCalc = cms.PSet(
    Class = cms.string("PileupGrid2d"),
    Grid2d = fftjet_pileup_grid_pf_calibrated,
    rhoFactor = cms.double(0.0)
)
fftjetPFJetMaker.pileupTableRecord = cms.string(fftjet_lut_types[etaSequenceTag].LUTRecord)
fftjetPFJetMaker.pileupTableName = cms.string("FFTPileupRhoEtaDependenceTable")
fftjetPFJetMaker.pileupTableCategory = cms.string("Pileup")
fftjetPFJetMaker.loadPileupFromDB = cms.bool(True)
# Configure module sequence
process.fftjetPFPatReco  = fftjetPFPatReco
process.fftjetPFJetReco  = fftjetPFJetMaker
#############end fftjet config################################
#############vplusganalyzer config############################
#make tags for LA modules
tags_la2 = cms.VInputTag()
tags_la3 = cms.VInputTag()
tags_la4 = cms.VInputTag()
rvalues  = cms.vdouble()
for R in calib_cone_radii:
					rtag  = "R%d" % (R*100.0)
					tags_la2.append("fftPFJetCorrectedLA2"+rtag)
					tags_la3.append("fftPFJetCorrectedLA3"+rtag)
					tags_la4.append("fftPFJetCorrectedLA4"+rtag)
					rvalues.append(R)


from TopQuarkAnalysis.VplusGAnalyzer.vplusganalyzer_fftjettest_cfi import * 
vplusganalyzer_fftjettest.srcFFTJetLA2 = tags_la2
vplusganalyzer_fftjettest.srcFFTJetLA3 = tags_la3
vplusganalyzer_fftjettest.srcFFTJetLA4 = tags_la4
vplusganalyzer_fftjettest.rValues      = rvalues

vplusganalyzer_fftjettest.srcPhoton    = cms.InputTag('photons')
vplusganalyzer_fftjettest.jetsForRho   = cms.string('kt6PFJets')
vplusganalyzer_fftjettest.runOnMC          = cms.bool(runOnMC)
vplusganalyzer_fftjettest.treeLabel          = cms.InputTag("fftjetPFPatReco", "FFTJetPatternRecognition")
vplusganalyzer_fftjettest.completeEventScale = cms.double(0.0)
vplusganalyzer_fftjettest.InitialScales      = fftjet_patreco_scales_50
vplusganalyzer_fftjettest.IsoValPhoton       = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                                                  cms.InputTag('phoPFIso:phIsoForGsfEle'),
                                                  cms.InputTag('phoPFIso:nhIsoForGsfEle'), )
vplusganalyzer_fftjettest.jetsForHT          = cms.InputTag('prunedGenParticles')
process.vplusganalyzer_fftjettest            = vplusganalyzer_fftjettest
#############vplusganalyzer config############################
process.options.wantSummary = cms.untracked.bool(wantSummary)
process.p = cms.Path(
                      process.makeCandidates*
											process.fftjetPFPatReco*
											process.fftjetPFJetReco
										)
######## make modules######################################## 
for la in la_jets:
	for R in calib_cone_radii:
			for config, modulename in adjust_cone_radius(R,fftjetPFJetMaker,la):
					setattr(process, modulename, config)
					process.p *= getattr(process, modulename)

process.p *= process.vplusganalyzer_fftjettest
#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('vplusg.root'),
#    outputCommands = cms.untracked.vstring('keep *')
#)
#process.e = cms.EndPath(process.out)

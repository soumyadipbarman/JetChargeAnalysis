import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

#process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/mc/RunIISummer20UL16MiniAODAPVv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/006C3FCB-7495-AF4F-A95F-B068A40DF498.root',
'/store/mc/RunIISummer20UL16MiniAODAPVv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/00B41AD6-70AC-3D47-807D-40F2D3744EC4.root',
#'/store/mc/RunIISummer20UL16MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/130000/1092DDD2-A0E6-DE47-8E72-807987BC5F99.root',
#'/store/mc/RunIISummer20UL16MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/140000/F2D918B0-03C3-EE4D-94AA-9EAA7FE84F19.root',
)

#eventsToSkip = cms.untracked.VEventRange('1:1950-1:2000'),
#eventsToSkip = cms.untracked.EventRange('1:1950-1:2000'),
#eventRanges = cms.untracked.VEventRange('1:1000-1:2000'),
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30000) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mcRun2_asymptotic_preVFP_v11')   #MC 2016APV
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mcRun2_asymptotic_v17')          #MC 2016

from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("HLTrigger.HLTcore.hltPrescaleRecorder_cfi")

process.options = cms.untracked.PSet(
)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)

process.MessageLogger = cms.Service("MessageLogger",  
 cout = cms.untracked.PSet(  
  default = cms.untracked.PSet( ## kill all messages in the log  
 
   limit = cms.untracked.int32(0)  
  ),  
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted  
 
   limit = cms.untracked.int32(-1)  
  )  
 ),  
 categories = cms.untracked.vstring('FwkJob'), 
 destinations = cms.untracked.vstring('cout')  
)  

# For Pileup JetID 
'''
process.load('RecoJets.JetProducers.PileupJetID_cfi')
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL17
process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL17),
    )
'''
process.load('RecoJets.JetProducers.PileupJetID_cfi')
#from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16APV     # APV
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16
process.pileupJetIdUpdated = process.pileupJetId.clone( 
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
	#algos = cms.VPSet(_chsalgos_106X_UL16APV)			     # APV
        algos = cms.VPSet(_chsalgos_106X_UL16)
    )

#For L1Prefiring
'''
from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
TheJets = cms.InputTag("slimmedJets"), #this should be the slimmedJets collection with up to date JECs !
DataEraECAL = cms.string("UL2016BtoH"),
DataEraMuon = cms.string("2016"),
UseJetEMPt = cms.bool(False),
PrefiringRateSystematicUnctyECAL = cms.double(0.2),
PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)
'''

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
TheJets = cms.InputTag("slimmedJets"), #this should be the slimmedJets collection with up to date JECs !
#DataEraECAL = cms.string("UL2016preVFP"),	# APV
#DataEraMuon = cms.string("2016preVFP"),
DataEraECAL = cms.string("UL2016postVFP"),
DataEraMuon = cms.string("2016postVFP"),
UseJetEMPt = cms.bool(False),
PrefiringRateSystematicUnctyECAL = cms.double(0.2),
PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)


#ak5 PF & Gen Jets
#from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
#from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
#from RecoMET.METProducers.PFMET_cfi import pfMet

#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')

# Select candidates that would pass CHS requirements
#process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#makes chs ak5 jets   (instead of ak4 that are default in miniAOD )
#process.ak5PFJetsCHS = ak5PFJets.clone(src = 'chs')

process.TFileService=cms.Service("TFileService",
    fileName=cms.string("Test_MC_2016UL.root")
)

print "test1"
# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
# Fix POWHEG if buggy (this PDF set will also appear on output,
# so only two more PDF sets can be added in PdfSetNames if not "")
#FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
#GenTag = cms.untracked.InputTag("genParticles"),
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring(
#               "CT10nlo.LHgrid"
                 "PDF4LHC15_100.LHgrid"
                , "MSTW2008nlo68cl.LHgrid"
                , "NNPDF30_100.LHgrid"
            )
      )

process.analyzeBasicPat = cms.EDAnalyzer("MiniAODAnalyzer",
#       photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
#       electronSrc = cms.untracked.InputTag("cleanPatElectrons"),
#       muonSrc = cms.untracked.InputTag("cleanPatMuons"),
#       tauSrc = cms.untracked.InputTag("cleanPatTaus"),
#        r = cms.EventRange('1:1000-1:2000'),
        jetSrc = cms.InputTag("slimmedJets"),
        metSrc = cms.InputTag("slimmedMETs"),
        genSrc = cms.untracked.InputTag("packedGenParticles"),
        pfSrc = cms.InputTag("packedPFCandidates"),
	metPATSrc = cms.InputTag("TriggerResults","","PAT"),
	metRECOSrc = cms.InputTag("TriggerResults","","RECO"),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
	objects = cms.InputTag("slimmedPatTrigger"),
        #objects = cms.InputTag("selectedPatTrigger"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	genparticles = cms.InputTag("prunedGenParticles"),
        bsSrc = cms.InputTag("offlineBeamSpot"),
        genjetSrc = cms.InputTag("slimmedGenJets"),
	genJetFlavourInfos = cms.InputTag('slimmedGenJetsFlavourInfos'),
        pileupSrc =cms.InputTag("slimmedAddPileupInfo"),
        ak5pfJetSrc = cms.InputTag("ak5PFJets"),
        ak5genJetSrc = cms.InputTag("ak5GenJets"),
        evtinfo =cms.InputTag("generator"),
        rho = cms.InputTag('fixedGridRhoAll'),
        LHEEventProductInputTag   = cms.InputTag('externalLHEProducer'),
        LHERunInfoProductInputTag = cms.InputTag('externalLHEProducer'),
        PDFCTEQWeightsInputTag   = cms.InputTag('pdfWeights:CT14'),
        PDFMMTHWeightsInputTag   = cms.InputTag('pdfWeights:MMHT2014lo68cl'),
        PDFNNPDFWeightsInputTag   = cms.InputTag('pdfWeights:NNPDF30'),
        #ak5PFJetCHSSrc    = cms.InputTag("ak5PFJetsCHS")
	RootFileName = cms.untracked.string('pythia8_test_13tev.root'),
	GenJET =  cms.untracked.bool(True),
	HistFill = cms.untracked.bool(True),
	MonteCarlo =  cms.untracked.bool(True),
	ParticleLabel =  cms.untracked.bool(False),
	Reconstruct =cms.untracked.bool(True),
#  EtaRange =  cms.untracked.double(5.0),
#  PtThreshold = cms.untracked.double(12.0),
  	EtaRange =  cms.untracked.double(3.0),
  	PtThreshold = cms.untracked.double(55.0), #effective is 21
  	LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81     
        tracks = cms.InputTag("generalTracks"), 
        filterGoodVertices = cms.InputTag("Flag_goodVertices"),
	filterglobalSuperTightHalo2016Filter  = cms.InputTag("Flag_globalSuperTightHalo2016Filter"),
	filterHBHENoiseFilter = cms.InputTag("Flag_HBHENoiseFilter"),
	filterHBHENoiseIsoFilter = cms.InputTag("Flag_HBHENoiseIsoFilter"),
	filterEcalDeadCellTriggerPrimitiveFilter = cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"),
	filterBadPFMuonFilter = cms.InputTag("Flag_BadPFMuonFilter"),
	filterBadPFMuonDzFilter = cms.InputTag("Flag_BadPFMuonDzFilter"),
	filtereeBadScFilter = cms.InputTag("Flag_eeBadScFilter"),
	filterecalBadCalibFilter = cms.InputTag("Flag_ecalBadCalibFilter"), 
#        scaleFactorsFile = cms.FileInPath('xxCondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('xxCondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#        scaleFactorsFile = cms.FileInPath('Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#        resolutionsFile = cms.FileInPath('Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
	#BTagEffFile = cms.string("BTagEfficiency2016preVFP_09Jun2024.root"),	# APV
	#BtagScaleFacFile = cms.string("btagging_2016preVFP.json.gz"),
	BTagEffFile = cms.string("BTagEfficiency2016postVFP_09Jun2024.root"),  # Non APV
	BtagScaleFacFile = cms.string("btagging_2016postVFP.json.gz"),
	bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
         'pfDeepCSVJetTags:probb',
         'pfDeepCSVJetTags:probbb',
	 'pfDeepFlavourJetTags:probb',
	 'pfDeepFlavourJetTags:probbb',
	 'pfDeepFlavourJetTags:problepb'
    )
 )


#process.options = cms.untracked.PSet(
#SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print "test2"
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)
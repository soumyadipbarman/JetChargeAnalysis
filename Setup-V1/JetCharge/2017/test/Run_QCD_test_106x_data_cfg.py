import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

#process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


# source
process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
'/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/260000/00EA7817-B2DB-E44A-9FD7-3FFB3778C024.root',
#'/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B6F8B6-90F1-E311-B72C-0025905A6092.root',
#'/store/data/Run2016F/JetHT/MINIAOD/21Feb2020_UL2016-v1/20000/0F71E50E-4F98-1446-B14A-A8600DDA6537.root',
#'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FFEEECDB-9680-6749-976A-F2CF385F8489.root',
#'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FEA84480-5EB1-DB42-B4C6-2F9229725E3C.root',
#'/store/data/Run2017C/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/00A34209-0E6A-9241-A53F-AB5EBCEB10F0.root',
#'/store/data/Run2017D/JetHT/MINIAOD/09Aug2019_UL2017-v1/00000/4EAE67FC-658E-1647-B91F-A978348DA634.root',
#'/store/data/Run2017E/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/016BE62D-1105-AC4F-8A58-59BD14326D8B.root',
#'/store/data/Run2017E/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/02CD90CE-110C-1B4C-BC0C-4BA10C64EFE2.root',
#'/store/data/Run2017F/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FF1B7616-113A-F14B-A97A-4939F16CB1B5.root',
#'/store/data/Run2017F/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FE81D290-531D-9044-BEC0-428C277EA7C6.root',
#'/store/data/Run2017F/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/FD6F6891-A1F9-534E-9B4C-6D49F8620D0B.root',
#'/store/data/Run2017B/JetHT/MINIAOD/31Mar2018-v1/90000/E87B6957-7139-E811-86AC-0CC47A78A468.root',
#'/store/data/Run2017B/JetHT/MINIAOD/31Mar2018-v1/90000/E6944CFD-8039-E811-BC15-0025905A606A.root',
#'/store/data/Run2017B/JetHT/MINIAOD/31Mar2018-v1/90000/E42325FE-8839-E811-A0EC-0CC47A78A436.root',
#'/store/data/Run2017C/JetHT/MINIAOD/31Mar2018-v1/30000/F019909F-4839-E811-A3BA-D4856445E5E2.root',
#'/store/data/Run2017E/JetHT/MINIAOD/31Mar2018-v1/90001/F41B041E-9837-E811-BFD2-008CFAC941DC.root',
#'/store/data/Run2017D/JetHT/MINIAOD/31Mar2018-v1/90000/E4402B21-F739-E811-BA23-0CC47A4C8E1E.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/02EDD8BA-B037-E811-B399-001E67792768.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/02FDDB8F-F037-E811-98EF-001E67792884.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/02FFDBF4-C037-E811-85A9-001E67E6F4EF.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8C3F3069-C937-E811-912E-A4BF0112DCB4.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8C0E44D5-1038-E811-9152-A4BF0112BE34.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8A922D57-C937-E811-9FFA-001E677926E4.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8A582EE4-CA37-E811-B066-A4BF0112BD18.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8A28CDED-ED37-E811-8584-002590796302.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8A05FC9C-AB37-E811-A5E2-A4BF01159734.root',
#'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/8A01C041-F137-E811-B683-001E677927F8.root',
#'/store/data/Run2017F/JetHT/MINIAOD/17Nov2017-v1/70000/F8D1C550-C9DF-E711-B91A-02163E0134D5.root',
#'/store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/159/00000/50A3A073-3B6C-E511-A997-02163E0144CD.root'
#'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/7E46D250-67B0-E511-BB96-0025905C3E66.root'
 )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.GlobalTag.globaltag = cms.string('POSTLS170_V5')
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_P_V56::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'GR_R_44_V11::All')
#process.GlobalTag = GlobalTag(process.GlobalTag,'74X_dataRun2_Prompt_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v2')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_ReReco_EOY17_v6')
#process.GlobalTag = GlobalTag(process.GlobalTag,'94X_dataRun2_v11')
#process.GlobalTag = GlobalTag(process.GlobalTag,'102X_dataRun2_v12')
#process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v28')
#process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v32')
#process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v35')
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_dataRun2_v37')

from PhysicsTools.PatAlgos.tools.coreTools import *
# produce PAT Layer 1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

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
process.load('RecoJets.JetProducers.PileupJetID_cfi')
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL17
process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL17),
    )

#process.load("HLTrigger.HLTcore.hltPrescaleRecorder_cfi")

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
    fileName=cms.string("Test_Data_2017UL.root")
)
print ("test1")
process.analyzeBasicPat = cms.EDAnalyzer("QCDEventShape",
#       photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
#       electronSrc = cms.untracked.InputTag("cleanPatElectrons"),
#       muonSrc = cms.untracked.InputTag("cleanPatMuons"),
#       tauSrc = cms.untracked.InputTag("cleanPatTaus"),
        jetSrc = cms.InputTag("slimmedJets"),
        metSrc = cms.InputTag("slimmedMETs"),
        genSrc = cms.untracked.InputTag("packedGenParticles"),
        pfSrc = cms.InputTag("packedPFCandidates"),
	metPATSrc = cms.InputTag("TriggerResults","","PAT"),
        metRECOSrc = cms.InputTag("TriggerResults","","RECO"),
        bits = cms.InputTag("TriggerResults","","HLT"),
        prescales = cms.InputTag("patTrigger"),
        objects = cms.InputTag("selectedPatTrigger"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	#genparticles = cms.InputTag("prunedGenParticles"),
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
	GenJET =  cms.untracked.bool(False),
	HistFill = cms.untracked.bool(True),
	MonteCarlo =  cms.untracked.bool(False),
	ParticleLabel =  cms.untracked.bool(False),
	Reconstruct =cms.untracked.bool(True),
#       EtaRange =  cms.untracked.double(5.0),
#       PtThreshold = cms.untracked.double(12.0),
  	EtaRange =  cms.untracked.double(3.0),
  	PtThreshold = cms.untracked.double(55.0), #effective is 21
  	LeadingPtThreshold = cms.untracked.double(150.0), #effective is 81       
#       scaleFactorsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'),
#       resolutionsFile = cms.FileInPath('CondFormats/JetMETObjects/data/Summer15_V0_MC_JER_AK4PFchs.txt'), 
#       scaleFactorsFile = cms.FileInPath('Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#       resolutionsFile = cms.FileInPath('Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
#       scaleFactorsFile = cms.FileInPath('Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
#       resolutionsFile = cms.FileInPath('Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
 	BTagEffFile = cms.string("BTagEfficiency2017_09Jun2024.root"),
	BtagScaleFacFile = cms.string("btagging_2017.json.gz"),
	bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
	 'pfDeepCSVJetTags:probb',
         'pfDeepCSVJetTags:probbb',
         'pfDeepFlavourJetTags:probb',
         'pfDeepFlavourJetTags:probbb',
         'pfDeepFlavourJetTags:problepb'
    )    
 )


#process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')
#process.analyzeBasicPat.append("keep *_ak5PFJets_*_EX")

#process.analyzeBasicPat.append("keep *_ak5PFJetsCHS_*_EX")
process.p = cms.Path(process.analyzeBasicPat)
print ("test2")
#process.p = cms.Path(process.ak5PFJets*process.ak5GenJets*process.analyzeBasicPat)

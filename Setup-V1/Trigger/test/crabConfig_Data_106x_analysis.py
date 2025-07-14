from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config # getUsernameFromSiteDB

config = Configuration()

config.section_("General")

                     ####### 22 DECEMBER 2020 ########

#config.General.requestName = 'Jet_Charge_JetHT_Run2016B-21Feb2020_ver1_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016B-21Feb2020_ver2_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016C-21Feb2020_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016D-21Feb2020_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016E-21Feb2020_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016F-21Feb2020_UL2016_HIPM-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016F-21Feb2020_UL2016-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016G-21Feb2020_UL2016-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016H-21Feb2020_UL2016-v1_22Dec'

#config.General.requestName = 'Jet_Charge_JetHT_Run2017B-09Aug2019_UL2017-v1_21Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017C-09Aug2019_UL2017-v1_21Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017D-09Aug2019_UL2017-v1_21Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017E-09Aug2019_UL2017-v1_21Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017F-09Aug2019_UL2017-v1_21Dec'

#config.General.requestName = 'Jet_Charge_JetHT_Run2018A-12Nov2019_UL2018-v2_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018B-12Nov2019_UL2018-v2_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018C-12Nov2019_UL2018_rsb-v1_22Dec'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018D-12Nov2019_UL2018_rsb-v1_22Dec'

                      ####### 12 JANUARY 2021 ########

#config.General.requestName = 'Jet_Charge_JetHT_Run2016B-21Feb2020_ver1_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016B-21Feb2020_ver2_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016C-21Feb2020_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016D-21Feb2020_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016E-21Feb2020_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016F-21Feb2020_UL2016_HIPM-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016F-21Feb2020_UL2016-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2016G-21Feb2020_UL2016-v1_12Jan'
config.General.requestName = 'Jet_Charge_JetHT_Run2016H-21Feb2020_UL2016-v1_12Jan'

#config.General.requestName = 'Jet_Charge_JetHT_Run2017B-09Aug2019_UL2017-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017C-09Aug2019_UL2017-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017D-09Aug2019_UL2017-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017E-09Aug2019_UL2017-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2017F-09Aug2019_UL2017-v1_12Jan'

#config.General.requestName = 'Jet_Charge_JetHT_Run2018A-12Nov2019_UL2018-v2_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018B-12Nov2019_UL2018-v2_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018C-12Nov2019_UL2018_rsb-v1_12Jan'
#config.General.requestName = 'Jet_Charge_JetHT_Run2018D-12Nov2019_UL2018_rsb-v1_12Jan'

config.General.workArea = 'crab_GJets'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'
#config.JobType.inputFiles= [
#"/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_MC_SF_AK4PFchs.txt", "/afs/cern.ch/work/t/tsarkar/private/QCD-13/CMSSW_7_6_3/src/Test/QCDEventShape/test/Fall15_25nsV2_DATA_UncertaintySources_AK4PF.txt"
#]

config.section_("Data")
#config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-21Feb2020_UL2016-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016G-21Feb2020_UL2016-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2016H-21Feb2020_UL2016-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017F-09Aug2019_UL2017-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2018A-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-12Nov2019_UL2018-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018C-12Nov2019_UL2018_rsb-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018D-12Nov2019_UL2018_rsb-v1/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.unitsPerJob = 200

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#config.Data.runRange = '246908-260627' # '193093-194075'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/'
config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'

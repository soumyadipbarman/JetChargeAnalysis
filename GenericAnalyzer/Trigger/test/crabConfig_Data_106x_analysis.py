from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config # getUsernameFromSiteDB

config = Configuration()
config.section_("General")

#config.General.requestName = 'JetHT_Run2016B-ver1_HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016B-ver2_HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016C-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016D-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016E-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016F-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016F_UL2016'
#config.General.requestName = 'JetHT_Run2016G_UL2016'
#config.General.requestName = 'JetHT_Run2016H_UL2016'

#config.General.requestName = 'JetHT_Run2017B_UL2017'
#config.General.requestName = 'JetHT_Run2017C_UL2017'
#config.General.requestName = 'JetHT_Run2017D_UL2017'
#config.General.requestName = 'JetHT_Run2017E_UL2017'
#config.General.requestName = 'JetHT_Run2017F_UL2017'

#config.General.requestName = 'JetHT_Run2018A_UL2018'
#config.General.requestName = 'JetHT_Run2018B_UL2018'
config.General.requestName = 'JetHT_Run2018C_UL2018'
#config.General.requestName = 'JetHT_Run2018D_UL2018'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/JetHT/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016G-UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016H-UL2016_MiniAODv2-v2/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018D-UL2018_MiniAODv2-v2/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
#config.Data.unitsPerJob = 200

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

#config.Data.runRange = '246908-260627' # '193093-194075'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

#config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/'
#config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/2016'
#config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/2017'
config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/2018'

#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016B-ver1_HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016B-ver2_HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016C-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016D-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016E-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016F-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016F_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016G_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016H_UL2016'

#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2017B_UL2017'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2017C_UL2017'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2017D_UL2017'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2017E_UL2017'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2017F_UL2017'

#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018A_UL2018'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018B_UL2018'
config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018C_UL2018'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018D_UL2018'

config.Data.publication = False
#config.Data.publishDataName = 'May2015_Data_analysis'

config.section_("Site")
#config.Site.storageSite = 'T2_IN_TIFR'
config.Site.storageSite = 'T3_US_FNALLPC'

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'JetHT_Run2016B-ver1_HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016B-ver2_HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016C-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016D-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016E-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016F-HIPM_UL2016'
#config.General.requestName = 'JetHT_Run2016F_UL2016'
#config.General.requestName = 'JetHT_Run2016G_UL2016'
#config.General.requestName = 'JetHT_Run2016H_UL2016'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'

config.JobType.inputFiles= [
"Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.txt",
"Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.txt",
"Summer19UL16APV_RunBCD_V7_DATA_UncertaintySources_AK4PFchs.txt",
"Summer19UL16APV_RunEF_V7_DATA_UncertaintySources_AK4PFchs.txt",
"Summer20UL16_JRV3_DATA_SF_AK4PFchs.txt",
"Summer20UL16_JRV3_DATA_PtResolution_AK4PFchs.txt",
"Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.txt"]

config.Data.inputDataset = '/JetHT/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016G-UL2016_MiniAODv2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016H-UL2016_MiniAODv2-v2/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'

config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Analysis/2016'
config.Data.publication = False

config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016B-ver1_HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016B-ver2_HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016C-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016D-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016E-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016F-HIPM_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016F_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016G_UL2016'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2016H_UL2016'

#config.Data.publishDataName = 'May2015_Data_analysis'

#config.Site.storageSite = 'T2_IN_TIFR'
config.Site.storageSite = 'T3_US_FNALLPC'

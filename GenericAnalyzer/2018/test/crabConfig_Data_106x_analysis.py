from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'JetHT_Run2018A_UL2018'
#config.General.requestName = 'JetHT_Run2018B_UL2018'
#config.General.requestName = 'JetHT_Run2018C_UL2018'
config.General.requestName = 'JetHT_Run2018D_UL2018'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'

config.JobType.inputFiles= ["Summer19UL18_JRV2_DATA_SF_AK4PFchs.txt", "Summer19UL18_JRV2_DATA_PtResolution_AK4PFchs.txt","Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.txt","Summer19UL18_RunB_V5_DATA_UncertaintySources_AK4PFchs.txt","Summer19UL18_RunC_V5_DATA_UncertaintySources_AK4PFchs.txt","Summer19UL18_RunD_V5_DATA_UncertaintySources_AK4PFchs.txt","BTagEfficiency2018_09Jun2024.root","btagging_2018.json.gz"]

#config.Data.inputDataset = '/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2018D-UL2018_MiniAODv2-v2/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Analysis/2018'
config.Data.publication = False

#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018A_UL2018'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018B_UL2018'
#config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018C_UL2018'
config.Data.outputDatasetTag = 'crab_projects_JetHT_Run2018D_UL2018'

#config.Data.publishDataName = 'May2015_Data_analysis'

#config.Site.storageSite = 'T2_IN_TIFR'
config.Site.storageSite = 'T3_US_FNALLPC'

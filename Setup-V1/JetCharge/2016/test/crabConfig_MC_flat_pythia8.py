#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

config.General.requestName ='JetCharge_TuneCP5_Flat_13TeV_pythia8_02072022'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run_QCD_test_106x_mc_cfg.py'

#config.JobType.psetName = options.cfg
#config.JobType.maxMemoryMB = 9000 # Default is 2500 : Max I have used is 13000
#config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000
#config.JobType.numCores = 4

config.JobType.inputFiles= ["/JER/MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt","/JER/MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt","/JER/MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt","/JER/MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt","/JEC/MC/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt","/JEC/MC/Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt"]

config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6_ext2-v2/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'EventBased'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 10  # for Automatic must be 180-2700 range
config.Data.unitsPerJob = 1  #For Filebased or Lumibased
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/%s/' % (sukundu)
config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Analysis'
config.Data.publication = False
config.Data.outputDatasetTag = 'JetCharge_MC_UL2017_Pythia8_flat_02072022'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite ='T2_IN_TIFR'
#config.Site.storageSite ='T3_US_FNALLPC'

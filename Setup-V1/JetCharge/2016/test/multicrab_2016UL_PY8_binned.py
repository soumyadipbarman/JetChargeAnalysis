#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

def getOptions():
        """
        Parse and return the arguments provided by the user.
        """
        usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
                 "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
                 "\nUse multicrab -h for help")
        parser = OptionParser(usage=usage)

        parser.add_option('-c', '--crabCmd',
                dest = 'crabCmd',
                default = '',
                help = "crab command",
                metavar = 'CMD')

        parser.add_option('-w', '--workArea',
                dest = 'workArea',
                default = '',
                help = "work area directory (only if CMD != 'submit')",
                metavar = 'WAD')

        parser.add_option('-m', '--crabCmdOpts',
                dest = 'crabCmdOpts',
                default = '',
                help = "options for crab command CMD",
                metavar = 'OPTS')

        (options, arguments) = parser.parse_args()



        if arguments:
                parser.error("Found positional argument(s): %s." % (arguments))
        if not options.crabCmd:
                parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
        if options.crabCmd != 'submit':
                if not options.workArea:
                        parser.error("(-w WAR, --workArea=WAR) option not provided.")
                if not os.path.isdir(options.workArea):
                        parser.error("'%s' is not a valid directory." % (options.workArea))
        return options

def main():
        options = getOptions()
# The submit command needs special treatment.
        if options.crabCmd == 'submit':
        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
                from CRABClient.UserUtilities import config
                config = config()

                config.General.requestName = None
                config.General.workArea = 'crab_projects'
                config.General.transferOutputs = True
                config.General.transferLogs = True
                config.JobType.pluginName = 'Analysis'
                config.JobType.psetName = 'Run_QCD_test_106x_mc_cfg.py'
	        #config.Data.ignoreLocality = True
                config.Data.inputDBS = 'global'
		config.JobType.inputFiles= ["Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt","Summer20UL16_JRV3_MC_SF_AK4PFchs.txt","Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt","btagging_2016postVFP.json.gz","BTagEfficiency2016postVFP_09Jun2024.root"]

		config.JobType.maxMemoryMB = 4000
		#config.JobType.priority = 9999
		#config.JobType.maxJobRuntimeMin = 3000
		config.JobType.numCores = 2
		#config.Data.splitting = 'EventAwareLumiBased'
                config.Data.splitting = 'FileBased'
                config.Data.unitsPerJob = 1
		#config.Data.useParent = True
                config.Data.inputDataset = None
		#config.Data.splitting = 'LumiBased'
		#config.Data.splitting = 'Automatic'
		#config.Data.unitsPerJob = 20
		#config.Data.totalUnits = 30
                config.Data.outputDatasetTag = None
		#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
		#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
		config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Analysis/2016'
                config.Data.publication = False
                config.JobType.allowUndistributedCMSSW = True
                #config.Site.storageSite = 'T2_IN_TIFR'
		config.Site.storageSite = 'T3_US_FNALLPC'

            # Will submit one task for each of these input datasets.
                inputDatasets = [
	'/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM',
	'/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM']
	
    	  	requestName =['Pt_15to30_TuneCP5','Pt_30to50_TuneCP5','Pt_50to80_TuneCP5','Pt_80to120_TuneCP5','Pt_120to170_TuneCP5','Pt_170to300_TuneCP5',
                               'Pt_300to470_TuneCP5','Pt_470to600_TuneCP5','Pt_600to800_TuneCP5','Pt_800to1000_TuneCP5','Pt_1000to1400_TuneCP5',
				'Pt_1400to1800_TuneCP5','Pt_1800to2400_TuneCP5','Pt_2400to3200_TuneCP5','Pt_3200toinf_TuneCP5']
		
                ireq = 0
                for inDS in inputDatasets: 
                        config.General.requestName = requestName[ireq]
                        ireq += 1
                        config.Data.inputDataset = inDS
                        config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
                        # Submit.
                        try:
                                print "Submitting for input dataset %s" % (inDS)
                                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
                        except HTTPException as hte:
                                print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
                        except ClientException as cle:
                                print "Submission for input dataset %s failed: %s" % (inDS, cle)
        elif options.workArea:
                for dir in os.listdir(options.workArea):
                        projDir = os.path.join(options.workArea, dir)
                        if not os.path.isdir(projDir):
                                continue
                        # Execute the crab command.
                        msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
                        print "-"*len(msg)
                        print msg
                        print "-"*len(msg)
                        try:
                                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
                        except HTTPException as hte:
                                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
                        except ClientException as cle:
                                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)
if __name__ == '__main__':
        main()

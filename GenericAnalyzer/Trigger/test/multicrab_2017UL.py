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
                config.JobType.allowUndistributedCMSSW = True
                config.JobType.psetName = 'Run_QCD_test_106x_data_cfg.py'
		#config.JobType.inputFiles= ["Summer19UL18_JRV2_MC_SF_AK4PFchs.txt","Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt","Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt","BTagEfficiency2018_09Jun2024.root","btagging_2018.json.gz"]

		#config.JobType.maxMemoryMB = 5000
		#config.JobType.priority = 9999
		#config.JobType.maxJobRuntimeMin = 3000
		#config.JobType.numCores = 2
		#config.Data.splitting = 'EventAwareLumiBased'
                #config.Data.splitting = 'FileBased'
                config.Data.inputDataset = None
                config.Data.inputDBS = 'global'
                #config.Data.ignoreLocality = True
                config.Data.splitting = 'LumiBased'
                config.Data.unitsPerJob = 15
		#config.Data.useParent = True
		#config.Data.totalUnits = 30
                config.Data.outputDatasetTag = None
		config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
		#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
		config.Data.outLFNDirBase = '/store/user/sobarman/JetCharge/Trigger/2017'
                config.Data.publication = False
                #config.Site.storageSite = 'T2_IN_TIFR'
		config.Site.storageSite = 'T3_US_FNALLPC'

            # Will submit one task for each of these input datasets.
                inputDatasets = [
		'/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD',
                '/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD',
                '/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD',
                '/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD',
		'/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD']

    	  	requestName =['JetHT_Run2017B_UL2017','JetHT_Run2017C_UL2017','JetHT_Run2017D_UL2017','JetHT_Run2017E_UL2017','JetHT_Run2017F_UL2017']
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

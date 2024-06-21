<< 'Test1'
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet60_v*" -o output_60.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet80_v*" -o output_80.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet140_v*" -o output_140.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet200_v*" -o output_200.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet260_v*" -o output_260.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet320_v*" -o output_320.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet400_v*" -o output_400.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet450_v*" -o output_450.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet500_v*" -o output_500.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet550_v*" -o output_550.csv
Test1

<< 'Test2'
cmsenv
pileupReCalc_HLTpaths.py -i output_60.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_60.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_80.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_80.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_140.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_140.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_200.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_200.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_260.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_260.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_320.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_320.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_400.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_400.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_450.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_450.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_500.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_500.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_550.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt -o HLT_PileupJSON_550.txt --runperiod Run2
Test2

#<< 'Test3'
cmsenv
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_60.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_60.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_80.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_80.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_140.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_140.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_200.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_200.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_260.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_260.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_320.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_320.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_400.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_400.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_450.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_450.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_500.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_500.root
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_550.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_550.root
#Test3

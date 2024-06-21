# JetChargeAnalysis

### How to run :
```
cmsrel  CMSSW_10_6_30
cd  CMSSW_10_6_30/src
cmsenv
mkdir Test
cd Test
mkedanlzr QCDEventShape
cd QCDEventShape
cd plugin
```
**In QCDEventShape folder, put the codes in plugin,test from those folders,
https://github.com/soumyadipbarman/JetChargeAnalysis/tree/master/JetCharge

 ### To complile and Run
 ```
 cd QCDEventShape
 scramv1 b    
 cd test 
 cmsRun Run_QCD_test_106x_mc_cfg.py    ## For MC
 cmsRun Run_QCD_test_106x_data_cfg.py  ## For DATA
```

### For TUnFold
 ```
1. Create a scram toolfile rootunfold.xml with content:
   <tool name="rootunfold" version="6.14.09">
   <info url="http://root.cern.ch/root/"/>
   <lib name="Unfold"/>
   <use name="roothistmatrix"/>
   <use name="rootxml"/>
   </tool>
2. Setup the new tool in your developer area
   scram setup /path/to/rootunfold.xml
3. Edit your Test/QCDEventShape/plugins/BuildFile.xml and add dependency on this new tool. All you need to add is the following line in your BuildFile.xml
  < use name="rootunfold"/>
4. Now rebuild: scram build

```

### Setup for MC Cross-Section Calculation
```
# activate McM tokens, must be done before setting cmsenv
cern-get-sso-cookie -u https://cms-pdmv.cern.ch/mcm/ -o ~/private/prod-cookie.txt --krb --reprocess
source /afs/cern.ch/cms/PPD/PdmV/tools/McM/getCookie.sh
# setup the grid xertificate
# grid-proxy-init -debug -verify
voms-proxy-init -voms cms
# setup cmssw release (choose the CMSSW version according to datasets)
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv

git cms-addpkg GeneratorInterface/Core
scram b -j8
cd ../../

cp -r calculateXSectionAndFilterEfficiency
cd calculateXSectionAndFilterEfficiency

# run using list of dataset names mode
./calculateXSectionAndFilterEfficiency.sh -f 2017/PY8Bin2017.txt -c RunIISummer20UL17MiniAODv2 -d MINIAODSIM -n 1000000

# For full details - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer
```

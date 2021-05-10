# JetChargeAnalysis

### How to run :
```
cmsrel  CMSSW_10_6_20
cd  CMSSW_10_6_20/src
cmsenv
mkdir Test
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
 cmsRun Run_QCD_test_106x_mc_cfg.py
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

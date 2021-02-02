# JetChargeAnalysis

### How to run :
```
cmsrel  CMSSW_10_2_18
cd  CMSSW_10_2_18/src
cmsenv
mkdir Test
mkedanlzr QCDEventShape
cd QCDEventShape
cd plugin
```
**In QCDEventShape folder, put the codes in plugin,test from those folders,
https://github.com/Sumankkundu/ChargedParticle/tree/master/QCDEventShape/2017/MC

 plugin : QCDEventShape.cc(Replace existing one) , EventShape_vector.cc(Copy)& EventShape_vector.h(copy) 
 
 test : All python files

 ### To complile and Run
 ```
 cd QCDEventShape
 scramv1 b    
 cd test 
 cmsRun Run_QCD_test_miaod_v2_76x_mc_cfg.py
```

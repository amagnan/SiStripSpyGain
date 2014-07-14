SiStripSpyGain
==============

#get release and package
ssh -X lxplus.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc472
cmsrel CMSSW_5_3_19
cd CMSSW_5_3_19/src
git cms-addpkg FWCore/Version
git cms-addpkg DQM/SiStripMonitorHardware
cd DQM/SiStripMonitorHardware/
scramv1 b -j 5

#get specific gain analysis
mkdir UserCode
cd UserCode
git clone 


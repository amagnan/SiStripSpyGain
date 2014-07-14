SiStripSpyGain
==============

#get release and package

ssh -X lxplus.cern.ch

export SCRAM_ARCH=slc6_amd64_gcc472

cmsrel CMSSW_5_3_19

cd CMSSW_5_3_19/src

git cms-addpkg FWCore/Version

git cms-addpkg DQM/SiStripMonitorHardware

scramv1 b -j 5

#get specific gain analysis

mkdir UserCode

cd UserCode

git clone https://github.com/amagnan/SiStripSpyGain.git

cd SiStripSpyGain/

scramv1 b -j 5

#prepare config, for example for spy run 198287 and global run 198372.

cd test/

#script prepareRunFolder.sh will create a local dir 198372 

# with config and file lists etc....

#source prepareRunFolder.sh <eos path to eos output area> \

# <spy run 3st 3 digits> <spy run last 3 digits> <global run> 

source prepareRunFolder.sh /store/group/comm_tracker/Strip/SpyChannel/ \

198 287 198372

cd 198372

#before running, files need to be staged or copied locally 

#(need enough space: copy locally only runs on which you will run often)

#to stage, execute script ./stageFiles_198372.sh.

#staging can take several hours

#to copy, execute script ./cpFiles_198372.sh

#copy is fast when files are already staged.

#run config interactively

cmsRun gainAnalysis_cfg.py

#to send in batch: edit and run ../submitBatchJob.sh


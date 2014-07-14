#!/bin/sh


cd /afs/cern.ch/work/a/amagnan/CMSSW_5_3_19/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh
RUNNUMBER=198372

eos mkdir -p /store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/
eos cp /afs/cern.ch/work/a/amagnan/CMSSW_5_3_19/src/UserCode/SiStripSpyGain/test/$RUNNUMBER/gainAnalysis_cfg.py /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/gainAnalysis_cfg.py

cp /afs/cern.ch/work/a/amagnan/CMSSW_5_3_19/src/UserCode/SiStripSpyGain/test/$RUNNUMBER/gainAnalysis_cfg.py .
cmsRun gainAnalysis_cfg.py

eos cp `pwd`/DQMStore.root /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/DQMStore.root
eos cp `pwd`/debug.log /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/debug.log
eos cp `pwd`/error.log /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/error.log
eos cp `pwd`/warning.log /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/warning.log
eos cp `pwd`/info.log /eos/cms/store/group/comm_tracker/Strip/SpyChannel/$RUNNUMBER/info.log

#echo "To submit batch job:"
#echo "bsub -q 1nd -J 198372 < submitBatchJob.sh"
#to request at least 30 GB space:
#bsub -R "pool>30000" -q 1nd -J job1 < runCastorJob.sh

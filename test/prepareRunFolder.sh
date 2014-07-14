#!/bin/sh


if (( "$#" != "4" ))
    then
    echo $# $*
    echo "Input parameter needed: <dest dir> <spy run first 3 digits> <spy run last 3 digits> <global run>"
    echo "castor path: /store/streamer/SiStripSpy/Commissioning11/000/198/287/"
    echo "eos path: /store/group/comm_tracker/Strip/SpyChannel/"
else
#INPUTFILE=$1
DESTDIR=$1
SPYRUNERA=$2
SPYRUNNUM=$3
RUNNUM=$4

source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh


INPUTDIR=root://castorcms//castor/cern.ch/cms/store/streamer/SiStripSpy/Commissioning11/000/$SPYRUNERA/$SPYRUNNUM/
echo "Input dir: " $INPUTDIR

mkdir -p $RUNNUM
eos ls $INPUTDIR/ | grep RUN00$RUNNUM > $RUNNUM/inputfiles.dat

echo "eos mkdir -p $DESTDIR/$RUNNUM/rawData/" > $RUNNUM/cpFiles_$RUNNUM.sh
awk '{$3="/eos/cms'$DESTDIR'/'$RUNNUM'/rawData/"$1}{$2="'$INPUTDIR'/"$1}{$1="eos cp "}{print $0}' $RUNNUM/inputfiles.dat >> $RUNNUM/cpFiles_$RUNNUM.sh

awk '{$1="stager_get -M '$INPUTDIR'/"$1}{print $0}' $RUNNUM/inputfiles.dat > $RUNNUM/stageFiles_$RUNNUM.sh

awk -v q="'" '{print q "'$DESTDIR'/'$RUNNUM'/rawData/"$1 q","}' $RUNNUM/inputfiles.dat > $RUNNUM/cfgFileList.txt

#cd /afs/cern.ch/work/a/amagnan/CMSSW_5_3_9/src/
#eval `scramv1 runtime -sh`
#cd -
#source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh

chmod u+x $RUNNUM/cpFiles_$RUNNUM.sh
chmod u+x $RUNNUM/stageFiles_$RUNNUM.sh
#./$RUNNUM/cpFiles_$RUNNUM.sh
#./$RUNNUM/stageFiles_$RUNNUM.sh

cat gainAnalysis_cfg_begin.py $RUNNUM/cfgFileList.txt gainAnalysis_cfg_end.py > $RUNNUM/gainAnalysis_cfg.py

#cat $RUNNUM/stageFiles_$RUNNUM.sh
#cat $RUNNUM/cpFiles_$RUNNUM.sh

fi
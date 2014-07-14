import FWCore.ParameterSet.Config as cms

process = cms.Process('SPYDQM')

#run 191226: 14/04/2012, 100 pb-1
#191247 78 pb-1
#200600: 09/08/2012 110 pb-1
#200525: 98 pb-1
#200473: 100 pb-1
#with spy...
#202209: 04/09/2012 17 pb-1 


#source of normal event data
process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(

    
       )
    )

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
  )

# --- Message Logging ---
#process.Tracer = cms.Service('Tracer',indentation = cms.untracked.string('$$'))
process.load('DQM.SiStripCommon.MessageLogger_cfi')
#process.MessageLogger.debugModules = cms.untracked.vstring('SiStripSpyMonitor')
process.MessageLogger.suppressInfo = cms.untracked.vstring('SiStripSpyDigiConverter')
process.MessageLogger.suppressWarning = cms.untracked.vstring('SiStripSpyDigiConverter')
#process.MessageLogger.suppressDebug = cms.untracked.vstring('SiStripSpyUnpacker')


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V43::All'

process.load('DQM.SiStripMonitorHardware.SiStripSpyUnpacker_cfi')
process.load('DQM.SiStripMonitorHardware.SiStripSpyDigiConverter_cfi')
process.SiStripSpyUnpacker.InputProductLabel = cms.InputTag('rawDataCollector')
process.SiStripSpyUnpacker.StoreScopeRawDigis = cms.bool(True)

process.SiStripSpyDigiConverter.InputProductLabel = cms.InputTag('SiStripSpyUnpacker','ScopeRawDigis')
process.SiStripSpyDigiConverter.StoreVirginRawDigis = cms.bool(True)
process.SiStripSpyDigiConverter.MinDigiRange = 100
process.SiStripSpyDigiConverter.MaxDigiRange = 1024
process.SiStripSpyDigiConverter.MinZeroLight = 0
process.SiStripSpyDigiConverter.MaxZeroLight = 1024
process.SiStripSpyDigiConverter.MinTickHeight = 0
process.SiStripSpyDigiConverter.MaxTickHeight = 1024
process.SiStripSpyDigiConverter.ExpectedPositionOfFirstHeaderBit = 6

# ---- FED Emulation ----
process.load('DQM.SiStripMonitorHardware.SiStripFEDEmulator_cfi')
process.SiStripFEDEmulator.SpyVirginRawDigisTag = cms.InputTag('SiStripSpyDigiConverter','VirginRaw')
process.SiStripFEDEmulator.ByModule = cms.bool(True) #use the digis stored by module (i.e. detId)

# ---- DQM
process.DQMStore = cms.Service("DQMStore")

process.load('UserCode.SiStripSpyGain.SiStripSpyGainAnalysis_cfi')
process.SiStripSpyGainAnalysis.SpyScopeRawDigisTag = cms.untracked.InputTag('SiStripSpyUnpacker','ScopeRawDigis')
process.SiStripSpyGainAnalysis.SpyPedSubtrDigisTag = cms.untracked.InputTag('SiStripFEDEmulator','PedSubtrModuleDigis')
process.SiStripSpyGainAnalysis.FillWithLocalEventNumber = False
process.SiStripSpyGainAnalysis.WriteDQMStore = True
process.SiStripSpyGainAnalysis.DQMStoreFileName = "DQMStore.root"
process.SiStripSpyGainAnalysis.SpyAPVeTag = cms.untracked.InputTag('SiStripSpyDigiConverter','APVAddress')
process.SiStripSpyGainAnalysis.ChannelList = cms.untracked.vuint32(470079220,369153672,369153636,369153668,369153640,369153676,369153660,369153644,369153656,402666258,402666261,402666257,402666262,369173524,369173540,369173556,369173560,369173544,369173528,369173532,369173548,369173564,369124374,369124381,369124373,369124378,369124377,369124382,436316292,436316296,436316300,436316304,436316308,436316312,436316316,436316320,436316324,436316328,436316332,436232454,436232466,436232469,436232474,436232457,436232462,436232465,436232470,436232473,436232453,436232458,436232461)

#TIB.1.4.1.1.1.* : 369173524,369173540,369173556,369173560,369173544,369173528,369173532,369173548,369173564
#TIB.1.1.1.1.1.* : 369124374,369124381,369124373,369124378,369124377,369124382
#TOB.1.6.1.1.1.* : 436316292,436316296,436316300,436316304,436316308,436316312,436316316,436316320,436316324,436316328,436316332
#TOB.1.1.1.1.1.* : 436232454,436232466,436232469,436232474,436232457,436232462,436232465,436232470,436232473,436232453,436232458,436232461

process.load('DQM.SiStripMonitorHardware.SiStripSpyMonitor_cfi')
process.SiStripSpyMonitor.SpyScopeRawDigisTag = cms.untracked.InputTag('SiStripSpyUnpacker','ScopeRawDigis')
process.SiStripSpyMonitor.SpyPedSubtrDigisTag = cms.untracked.InputTag('SiStripFEDEmulator','PedSubtrModuleDigis')
process.SiStripSpyMonitor.SpyAPVeTag = cms.untracked.InputTag('SiStripSpyDigiConverter','APVAddress')
process.SiStripSpyMonitor.FillWithLocalEventNumber = False
process.SiStripSpyMonitor.WriteDQMStore = True
process.SiStripSpyMonitor.DQMStoreFileName = "DQMStore.root"
#process.SiStripSpyMonitor.OutputErrors = "NoData","MinZero","MaxSat","LowRange","HighRange","LowDAC","HighDAC","OOS","OtherPbs","APVError","APVAddressError","NegPeds"
#process.SiStripSpyMonitor.OutputErrors = "MinZero","MaxSat","LowRange","HighRange","LowDAC","HighDAC","OOS","OtherPbs","APVError","APVAddressError","NegPeds"
#process.SiStripSpyMonitor.WriteCabling = True

#---- TkHistoMap -----
#process.load("DQM.SiStripCommon.TkHistoMap_cfi")
process.TkDetMap = cms.Service("TkDetMap")
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")

#process.load('PerfTools.Callgrind.callgrindSwitch_cff')

process.p = cms.Path(
    process.SiStripSpyUnpacker
    *process.SiStripSpyDigiConverter
    #*process.profilerStart
    *process.SiStripFEDEmulator
    *process.SiStripSpyGainAnalysis
    *process.SiStripSpyMonitor
    #*process.profilerStop 
    )

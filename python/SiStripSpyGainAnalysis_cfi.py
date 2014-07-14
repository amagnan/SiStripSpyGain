import FWCore.ParameterSet.Config as cms


SiStripSpyGainAnalysis = cms.EDAnalyzer(
    "SiStripSpyGainAnalysis",
    #Raw data collection
    SpyScopeRawDigisTag = cms.untracked.InputTag('SiStripSpyUnpacker','ScopeRawDigis'),
    SpyPedSubtrDigisTag = cms.untracked.InputTag('SiStripFEDEmulator','PedSubtrModuleDigis'),
    SpyL1Tag = cms.untracked.InputTag('SiStripSpyUnpacker','L1ACount'),
    SpyTotalEventCountTag = cms.untracked.InputTag('SiStripSpyUnpacker','TotalEventCount'),
    SpyAPVeTag = cms.untracked.InputTag('SiStripSpyDigiConverter','APVAddress'),
    MinDigiRange = cms.untracked.uint32(400),
    MaxDigiRange = cms.untracked.uint32(950),
    MinZeroLight = cms.untracked.uint32(15),
    MaxZeroLight = cms.untracked.uint32(180),
    MinTickHeight = cms.untracked.uint32(555),
    MaxTickHeight = cms.untracked.uint32(1010),
    #Folder in DQM Store to write global histograms to
    HistogramFolderName = cms.untracked.string('SiStrip/ReadoutView/SpyMonitoringSummary'),
    #Fill all detailed histograms at FED level even if they will be empty (so that files can be merged)
    FillAllDetailedHistograms = cms.untracked.bool(False),
    FillWithEventNumber = cms.untracked.bool(True),
    FillWithLocalEventNumber = cms.untracked.bool(False),
    #Whether to write the DQM store to a file at the end of the run and the file name
    WriteDQMStore = cms.untracked.bool(True),
    DQMStoreFileName = cms.untracked.string('DQMStore.root'),
    #OutputErrors = cms.untracked.vstring('NoData','MinZero','MaxSat','LowRange','HighRange','LowDAC','HighDAC','OOS','OtherPbs','APVError','APVAddressError','NegPeds'),
    OutputErrors = cms.untracked.vstring(),
    WriteCabling = cms.untracked.bool(False),
    ChannelList = cms.untracked.vuint32(),
    #Histogram configuration
    ZeroLightvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                         NBins = cms.untracked.uint32(600),
                                                         Min = cms.untracked.double(0),
                                                         Max = cms.untracked.double(600) ),
    
    TickHeightvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),

    RatiovsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),

    GainvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),

    DBGainvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),

    GainDiffvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),

    RelGainDiffvsTimeHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(600),
                                                             Min = cms.untracked.double(0),
                                                             Max = cms.untracked.double(600) ),
    
    allEvtAllChGainDiffHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(100),
                                                             Min = cms.untracked.double(-2),
                                                             Max = cms.untracked.double(2) ),

    allEvtAllChGainRatioHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                                             NBins = cms.untracked.uint32(100),
                                                             Min = cms.untracked.double(-2),
                                                             Max = cms.untracked.double(2) ),
    


    RunMeanGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),

    RunMeanDBGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),
    

    RunRmsGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),

    RunRmsDBGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),

    RunMeanGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),

    
    RunMeanRelGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),
    
    RunRmsGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),

    RunRmsRelGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True),
                                NBins = cms.untracked.uint32(100),
                                 Min = cms.untracked.double(-2),
                                 Max = cms.untracked.double(2) ),



    #2D Histos
    allEvtAllChGainScatterHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),

    fedIdvsChannelIdRunMeanZeroLightPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanTickHeightPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanRatioPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanDBGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunMeanRelGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsZeroLightPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsTickHeightPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsRatioPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsDBGainPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),
    fedIdvsChannelIdRunRmsRelGainDiffPerChHistogramConfig = cms.untracked.PSet( Enabled = cms.untracked.bool(True)),

     TkHistoMapHistogramConfig = cms.untracked.PSet( 
         Enabled = cms.untracked.bool(True) 
         ),

)

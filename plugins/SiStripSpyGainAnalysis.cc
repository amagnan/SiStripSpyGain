// Original Author:  Anne-Marie Magnan
//         Created:  2010/01/11
//         Update :  2013/05/06
//

#include <sstream>
#include <memory>
#include <list>
#include <algorithm>
#include <cassert>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForHardwareSystems.h"

#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
#include "CondFormats/DataRecord/interface/SiStripPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"

#include "EventFilter/SiStripRawToDigi/interface/SiStripFEDBuffer.h"

// For plotting.
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DQM/SiStripMonitorHardware/interface/SiStripFEDSpyBuffer.h"
#include "DQM/SiStripMonitorHardware/interface/SiStripSpyUtilities.h"
#include "UserCode/SiStripSpyGain/interface/GMHistograms.h"


// For DB

#include "CondFormats/DataRecord/interface/SiStripApvGainRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"



//
// Class declaration
//

class SiStripSpyGainAnalysis : public edm::EDAnalyzer
{
 public:


  explicit SiStripSpyGainAnalysis(const edm::ParameterSet&);
  ~SiStripSpyGainAnalysis();

 private:

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();


  //tag of spydata collection
  edm::InputTag spyScopeRawDigisTag_;

  //tag of l1A and apveAddress counters
  edm::InputTag spyL1Tag_;
  edm::InputTag spyTotCountTag_;
  edm::InputTag spyAPVeTag_;

  uint32_t minDigiRange_;
  uint32_t maxDigiRange_;
  uint32_t minDigitalLow_;
  uint32_t maxDigitalLow_;
  uint32_t minDigitalHigh_;
  uint32_t maxDigitalHigh_;

  unsigned int evt_;

  DQMStore* dqm_;
  //folder name for histograms in DQMStore
  std::string folderName_;
  //do histos vs time with time=event number. Default time = orbit number (s)
  bool fillWithEvtNum_;
  bool fillWithLocalEvtNum_;
  //write the DQMStore to a root file at the end of the job
  bool writeDQMStore_;
  std::string dqmStoreFileName_;

  GMHistograms histManager_;  
  uint16_t firstHeaderBit_;

  sistrip::SpyUtilities utility_;
  sistrip::SpyUtilities::FrameQuality frameQuality_;

  GMHistograms::Statistics channelData_[sistrip::FED_ID_MAX-sistrip::FED_ID_MIN+1][sistrip::FEDCH_PER_FED];

  bool isValid_[sistrip::FED_ID_MAX-sistrip::FED_ID_MIN+1][sistrip::FEDCH_PER_FED];
  std::map<uint32_t,GMHistograms::Statistics> channelMap_;

  const SiStripFedCabling* lCabling;

  std::vector<unsigned> channelList_;


};

using edm::LogError;
using edm::LogWarning;
using edm::LogInfo;
//
// Constructors and destructor
//

SiStripSpyGainAnalysis::SiStripSpyGainAnalysis(const edm::ParameterSet& iConfig)
  : spyScopeRawDigisTag_(iConfig.getUntrackedParameter<edm::InputTag>("SpyScopeRawDigisTag",edm::InputTag("SiStripSpyChannelScopeDigis","ScopeRawDigis"))),
    spyL1Tag_(iConfig.getUntrackedParameter<edm::InputTag>("SpyL1Tag",edm::InputTag("SiStripSpyChannelDataMerger","SpyL1ACount"))),
    spyTotCountTag_(iConfig.getUntrackedParameter<edm::InputTag>("SpyTotalEventCountTag",edm::InputTag("SiStripSpyChannelDataMerger","SpyTotalEventCount"))),
    spyAPVeTag_(iConfig.getUntrackedParameter<edm::InputTag>("SpyAPVeTag",edm::InputTag("SiStripSpyChannelDataMerger","SpyAPVAddress"))),
    dqm_(0),
    folderName_(iConfig.getUntrackedParameter<std::string>("HistogramFolderName","SiStrip/ReadoutView/SpyMonitoringSummary")),
    fillWithEvtNum_(iConfig.getUntrackedParameter<bool>("FillWithEventNumber",false)),
    fillWithLocalEvtNum_(iConfig.getUntrackedParameter<bool>("FillWithLocalEventNumber",false)),
    writeDQMStore_(iConfig.getUntrackedParameter<bool>("WriteDQMStore",false)),
    dqmStoreFileName_(iConfig.getUntrackedParameter<std::string>("DQMStoreFileName","DQMStore.root")),
    firstHeaderBit_(0),
    channelList_(iConfig.getUntrackedParameter<std::vector<unsigned> >("ChannelList"))
{

  evt_ = 0;
  std::ostringstream pDebugStream;
  histManager_.initialise(iConfig,&pDebugStream);

  frameQuality_.minDigiRange = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MinDigiRange",100));
  frameQuality_.maxDigiRange = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MaxDigiRange",1024));
  frameQuality_.minZeroLight = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MinZeroLight",0));
  frameQuality_.maxZeroLight = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MaxZeroLight",1024));
  frameQuality_.minTickHeight = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MinTickHeight",0));
  frameQuality_.maxTickHeight = static_cast<uint16_t>(iConfig.getUntrackedParameter<uint32_t>("MaxTickHeight",1024));
  
}

SiStripSpyGainAnalysis::~SiStripSpyGainAnalysis()
{ 

}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripSpyGainAnalysis::beginJob()
{

  //get DQM store
  dqm_ = &(*edm::Service<DQMStore>());
  dqm_->setCurrentFolder(folderName_);

  LogInfo("SiStripSpyGainAnalysis") << " Histograms will be written in " 
				     << folderName_ 
				     << ". Current folder is : " 
				     << dqm_->pwd() 
				     << std::endl;
  
  //this propagates dqm_ to the histoclass, must be called !
  histManager_.bookTopLevelHistograms(dqm_);
  
  evt_ = 0;
  firstHeaderBit_ = 0;


  for (uint16_t lFedId = sistrip::FED_ID_MIN; lFedId <= sistrip::FED_ID_MAX; ++lFedId) {//loop on feds
    for (uint16_t lFedChannel = 0; lFedChannel < sistrip::FEDCH_PER_FED; lFedChannel++){//loop on channels
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh = 0;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh = 0;


      isValid_[lFedId-sistrip::FED_ID_MIN][lFedChannel] = false;
    }
  }

}



// ------------ method called to for each event  ------------
void
SiStripSpyGainAnalysis::analyze(const edm::Event& iEvent, 
				  const edm::EventSetup& iSetup)
{

  static bool firstEvent = true;

   //update gains
  edm::ESHandle<SiStripApvGain> lgainHandle = utility_.getGainHandle(iSetup);

  //update cabling and pedestals
  lCabling = utility_.getCabling( iSetup );

  //For spy data
  //get map of TotalEventCount and L1ID, indexed by fedId, and APVaddress indexed by fedIndex.
  edm::Handle<std::vector<uint32_t> > lSpyL1IDHandle,lSpyTotCountHandle,lSpyAPVeHandle;
  try {
    iEvent.getByLabel(spyL1Tag_,lSpyL1IDHandle);
    iEvent.getByLabel(spyTotCountTag_,lSpyTotCountHandle);
    iEvent.getByLabel(spyAPVeTag_,lSpyAPVeHandle);
  }
  catch (const cms::Exception& e) {
    LogError("SiStripSpyGainAnalysis") << e.what() ;
    return;
  }
  //const std::vector<uint32_t> & lSpyMaxCountMap = *lSpyL1IDHandle;
  //const std::vector<uint32_t> & lSpyMinCountMap = *lSpyTotCountHandle;
  const std::vector<uint32_t> & lSpyAPVeVec = *lSpyAPVeHandle;

  //retrieve the scope digis
  edm::Handle<edm::DetSetVector<SiStripRawDigi> > digisHandle;
  try {
    iEvent.getByLabel(spyScopeRawDigisTag_, digisHandle);
  }
  catch (const cms::Exception& e) {
    LogError("SiStripSpyGainAnalysis") << e.what() ;
    return;
  }
  const edm::DetSetVector<SiStripRawDigi> * lInputDigis = digisHandle.product();


  //for first event only
  //loop first on channels to calculate majority value of when the first header bit is found.
  //output info message to give the value found
  //should currently be 6 but may vary in the futur
  //then we can check firstTrailerBit is +256+24 after

  if (firstEvent){
    sistrip::SpyUtilities::getMajorityHeader(lInputDigis,firstHeaderBit_);
  }


  double lTime;
  //if (fillWithEvtNum_) 
  //lTime = iEvent.id().event();
  //else if (fillWithLocalEvtNum_) lTime = evt_;
  //no orbit number for spy data !!
  //else lTime = iEvent.orbitNumber()/11223.;
  lTime = iEvent.id().event();
  if (fillWithLocalEvtNum_) lTime = evt_;


  //std::cout << " -- Processing event " << evt_ << std::endl;

  //loop over all FEDs and channels
  
  for (uint16_t lFedId = sistrip::FED_ID_MIN; lFedId <= sistrip::FED_ID_MAX; ++lFedId) {//loop on feds

    unsigned int lAPVAddrRef = lSpyAPVeVec.at(lFedId);

    //uint32_t lDetIdPrevious = -1;
    //size_t apvCounter=0;

    for (uint16_t lFedChannel = 0; lFedChannel < sistrip::FEDCH_PER_FED; lFedChannel++){//loop on channels
      
      uint32_t lFedIndex = sistrip::FEDCH_PER_FED*lFedId + lFedChannel;
      
      const FedChannelConnection & lConnection = lCabling->connection(lFedId,lFedChannel);

      if (!lConnection.isConnected()) continue;

      uint32_t lDetId = lConnection.detId();
      
      //uint16_t lNPairs = lConnection.nApvPairs();
      uint16_t lPair = lConnection.apvPairNumber();

      edm::DetSetVector<SiStripRawDigi>::const_iterator lDigis = lInputDigis->find(lFedIndex);
       
      //no digis found, continue.
      if (lDigis == lInputDigis->end()) {
	LogDebug("SiStripSpyGainAnalysis") << " -- digis not found in ScopeRawDigis map for FEDID " 
						  << lFedId << " and FED channel " << lFedChannel << std::endl;
	continue;
      }

      sistrip::SpyUtilities::Frame lFrame = sistrip::SpyUtilities::extractFrameInfo(*lDigis);
      

      if (!sistrip::SpyUtilities::isValid(lFrame,frameQuality_,firstHeaderBit_)) continue;
      if (lFrame.apvAddress.first != lAPVAddrRef || 
	  lFrame.apvAddress.second != lAPVAddrRef) continue;

      isValid_[lFedId-sistrip::FED_ID_MIN][lFedChannel] = true;


      // For DB
      SiStripApvGain::Range lgainRange = lgainHandle->getRange(lDetId);
      SiStripApvGain lapvGain;
      lapvGain.put(lDetId,lgainRange);
      SiStripGain gain(lapvGain,1.);
      //size_t apvE= (lgainRange.second-lgainRange.first);

      float lDBGain = gain.getApvGain(2*lPair,lgainRange);

      //if(lDetIdPrevious != lDetId){apvCounter=0;  lDBGain=gain.getApvGain(apvCounter,lgainRange);}
      //if(lDetIdPrevious == lDetId){apvCounter+=2; lDBGain=gain.getApvGain(apvCounter,lgainRange);}
      //lDetIdPrevious = lDetId;
 
      for (unsigned iCh(0); iCh<channelList_.size(); ++iCh){
	if (lDetId == channelList_[iCh]){
	  GMHistograms::Quantities lQuant;
	  lQuant.ZeroLight = lFrame.digitalLow;
	  lQuant.TickHeight = lFrame.digitalHigh;
	  lQuant.Ratio = lFrame.Ratio();
	  lQuant.Gain = lFrame.Gain();
	  lQuant.DBGain = lDBGain;
	  lQuant.GainDiff = lQuant.Gain - lDBGain ;
	  lQuant.RelGainDiff = (lQuant.Gain - lDBGain)/ lDBGain ;
	  histManager_.fillVsTimeHistograms(iCh,lPair,lQuant,lTime);
	  //std::cout << " -- Processing detid " << lDetId << std::endl;

	  break;
	}
      }

      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh += lFrame.digitalLow;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh += lFrame.digitalLow*lFrame.digitalLow;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh++;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh += lFrame.digitalHigh;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh += lFrame.digitalHigh*lFrame.digitalHigh;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh++;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh += lFrame.Ratio();
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh += lFrame.Ratio()*lFrame.Ratio();
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh++;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh += lFrame.Gain() ;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh += lFrame.Gain()*lFrame.Gain();
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh++;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh += lDBGain;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh += lDBGain * lDBGain;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh++;

      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh += lFrame.Gain()- lDBGain ;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh += ( lFrame.Gain() - lDBGain)
	                                                                    *( lFrame.Gain() - lDBGain );
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh++;


      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh += (lFrame.Gain()- lDBGain)/lDBGain;
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh += ( (lFrame.Gain() - lDBGain)/lDBGain)
                                                                                     *( (lFrame.Gain() - lDBGain)/lDBGain );
      channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh++;


      std::pair<std::map<uint32_t,GMHistograms::Statistics>::iterator,bool> lIsInserted = 
	channelMap_.insert(std::pair<uint32_t,GMHistograms::Statistics>(lDetId,channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel]));
      //if already exists, increment statistics counters
      if (!lIsInserted.second) {
	GMHistograms::Statistics & lStats = (lIsInserted.first)->second;
	lStats.RunMeanZeroLightPerCh += lFrame.digitalLow;
	lStats.RunRmsZeroLightPerCh += lFrame.digitalLow*lFrame.digitalLow;
	lStats.CounterZeroLightPerCh++;
	lStats.RunMeanTickHeightPerCh += lFrame.digitalHigh;
	lStats.RunRmsTickHeightPerCh += lFrame.digitalHigh*lFrame.digitalHigh;
	lStats.CounterTickHeightPerCh++;
	lStats.RunMeanRatioPerCh += lFrame.Ratio();
	lStats.RunRmsRatioPerCh += lFrame.Ratio()*
                                   lFrame.Ratio();
	lStats.CounterRatioPerCh++;
	lStats.RunMeanGainPerCh += lFrame.Gain();
	lStats.RunRmsGainPerCh += lFrame.Gain()*
                                   lFrame.Gain();
	lStats.CounterGainPerCh++;
	lStats.RunMeanDBGainPerCh += lDBGain;
	lStats.RunRmsDBGainPerCh += lDBGain * lDBGain;
	lStats.CounterDBGainPerCh++;
	lStats.RunMeanGainDiffPerCh += lFrame.Gain() - lDBGain;
	lStats.RunRmsGainDiffPerCh += ( lFrame.Gain() - lDBGain )*
                                   ( lFrame.Gain() - lDBGain);
	lStats.CounterGainDiffPerCh++;

	lStats.RunMeanRelGainDiffPerCh += (lFrame.Gain() - lDBGain)/lDBGain;
	lStats.RunRmsRelGainDiffPerCh += ( (lFrame.Gain() - lDBGain)/lDBGain )*
                                      ( (lFrame.Gain() - lDBGain)/lDBGain);
	lStats.CounterRelGainDiffPerCh++;
      }

      histManager_.fillAllEvtAllChHistograms(lDBGain, lFrame.Gain());


    }//loop on channels
  }//loop on feds

  //used to fill histo vs time with local event number....
  evt_++;
  firstEvent = false;
  
}//analyze method






// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripSpyGainAnalysis::endJob()
{

  //loop over all FEDs and channels
  
  for (uint16_t lFedId = sistrip::FED_ID_MIN; lFedId <= sistrip::FED_ID_MAX; ++lFedId) {//loop on feds

    for (uint16_t lFedChannel = 0; lFedChannel < sistrip::FEDCH_PER_FED; lFedChannel++){//loop on channels

      //const FedChannelConnection & lConnection = lCabling->connection(lFedId,lFedChannel);

      //if (!lConnection.isConnected()) continue;

      //uint32_t lDetId = lConnection.detId();

      if (!isValid_[lFedId-sistrip::FED_ID_MIN][lFedChannel]) continue;

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh 
             - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh)>=0){
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterZeroLightPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh = 0;
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanZeroLightPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsZeroLightPerCh = 0;
      }

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterTickHeightPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanTickHeightPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsTickHeightPerCh = 0;
      }

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRatioPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRatioPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRatioPerCh = 0;
      }


      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainPerCh = 0;
      }

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterDBGainPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanDBGainPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsDBGainPerCh = 0;
      }

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterGainDiffPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanGainDiffPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsGainDiffPerCh = 0;
      }

      if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh != 0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh =
            channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh;
         if(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh 
	  - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh) >=0){
         channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh =
          sqrt(  channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh / channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].CounterRelGainDiffPerCh 
	      - (channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh)*(channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh));
         }else{
          channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh = 0; 
         }
      }else{
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunMeanRelGainDiffPerCh = 0;
	 channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel].RunRmsRelGainDiffPerCh = 0;
      }

      histManager_.fillRunMeanAndRunRmsHistograms(lFedChannel,lFedId,channelData_[lFedId-sistrip::FED_ID_MIN][lFedChannel]);
    }
  }

  std::map<uint32_t,GMHistograms::Statistics>::iterator lMapIter = channelMap_.begin();
  for (; lMapIter != channelMap_.end(); lMapIter++){

    GMHistograms::Statistics lStat = lMapIter->second;
    uint32_t lDetId = lMapIter->first;

    if(lStat.CounterZeroLightPerCh !=0){
      lStat.RunMeanZeroLightPerCh = lStat.RunMeanZeroLightPerCh / lStat.CounterZeroLightPerCh;
      if(lStat.RunRmsZeroLightPerCh / lStat.CounterZeroLightPerCh - lStat.RunMeanZeroLightPerCh*lStat.RunMeanZeroLightPerCh){
       lStat.RunRmsZeroLightPerCh = sqrt( lStat.RunRmsZeroLightPerCh / lStat.CounterZeroLightPerCh - lStat.RunMeanZeroLightPerCh*lStat.RunMeanZeroLightPerCh );
      }else{
       lStat.RunRmsZeroLightPerCh = 0;
      }
    }else{
      lStat.RunMeanZeroLightPerCh = 0;
      lStat.RunRmsZeroLightPerCh = 0;
    }

    if(lStat.CounterTickHeightPerCh !=0){
      lStat.RunMeanTickHeightPerCh = lStat.RunMeanTickHeightPerCh / lStat.CounterTickHeightPerCh;
      if(lStat.RunRmsTickHeightPerCh / lStat.CounterTickHeightPerCh - lStat.RunMeanTickHeightPerCh*lStat.RunMeanTickHeightPerCh){
       lStat.RunRmsTickHeightPerCh = sqrt( lStat.RunRmsTickHeightPerCh / lStat.CounterTickHeightPerCh - lStat.RunMeanTickHeightPerCh*lStat.RunMeanTickHeightPerCh );
      }else{
       lStat.RunRmsTickHeightPerCh = 0;
      }
    }else{
      lStat.RunMeanTickHeightPerCh = 0;
      lStat.RunRmsTickHeightPerCh = 0;
    }

    if(lStat.CounterRatioPerCh !=0){
      lStat.RunMeanRatioPerCh = lStat.RunMeanRatioPerCh / lStat.CounterRatioPerCh;
      if(lStat.RunRmsRatioPerCh / lStat.CounterRatioPerCh - lStat.RunMeanRatioPerCh*lStat.RunMeanRatioPerCh){
       lStat.RunRmsRatioPerCh = sqrt( lStat.RunRmsRatioPerCh / lStat.CounterRatioPerCh - lStat.RunMeanRatioPerCh*lStat.RunMeanRatioPerCh );
      }else{
       lStat.RunRmsRatioPerCh = 0;
      }
    }else{
      lStat.RunMeanRatioPerCh = 0;
      lStat.RunRmsRatioPerCh = 0;
    }

    if(lStat.CounterGainPerCh !=0){
      lStat.RunMeanGainPerCh = lStat.RunMeanGainPerCh / lStat.CounterGainPerCh;
      if(lStat.RunRmsGainPerCh / lStat.CounterGainPerCh - lStat.RunMeanGainPerCh*lStat.RunMeanGainPerCh){
       lStat.RunRmsGainPerCh = sqrt( lStat.RunRmsGainPerCh / lStat.CounterGainPerCh - lStat.RunMeanGainPerCh*lStat.RunMeanGainPerCh );
      }else{
       lStat.RunRmsGainPerCh = 0;
      }
    }else{
      lStat.RunMeanGainPerCh = 0;
      lStat.RunRmsGainPerCh = 0;
    }


    if(lStat.CounterDBGainPerCh !=0){
      lStat.RunMeanDBGainPerCh = lStat.RunMeanDBGainPerCh / lStat.CounterDBGainPerCh;
      if(lStat.RunRmsDBGainPerCh / lStat.CounterDBGainPerCh - lStat.RunMeanDBGainPerCh*lStat.RunMeanDBGainPerCh){
       lStat.RunRmsDBGainPerCh = sqrt( lStat.RunRmsDBGainPerCh / lStat.CounterDBGainPerCh - lStat.RunMeanDBGainPerCh*lStat.RunMeanDBGainPerCh );
      }else{
       lStat.RunRmsDBGainPerCh = 0;
      }
    }else{
      lStat.RunMeanDBGainPerCh = 0;
      lStat.RunRmsDBGainPerCh = 0;
    }

    if(lStat.CounterGainDiffPerCh !=0){
      lStat.RunMeanGainDiffPerCh = lStat.RunMeanGainDiffPerCh / lStat.CounterGainDiffPerCh;
      if(lStat.RunRmsGainDiffPerCh / lStat.CounterGainDiffPerCh - lStat.RunMeanGainDiffPerCh*lStat.RunMeanGainDiffPerCh){
       lStat.RunRmsGainDiffPerCh = sqrt( lStat.RunRmsGainDiffPerCh / lStat.CounterGainDiffPerCh - lStat.RunMeanGainDiffPerCh*lStat.RunMeanGainDiffPerCh );
      }else{
       lStat.RunRmsGainDiffPerCh = 0;
      }
    }else{
      lStat.RunMeanGainDiffPerCh = 0;
      lStat.RunRmsGainDiffPerCh = 0;
    }

    if(lStat.CounterRelGainDiffPerCh !=0){
      lStat.RunMeanRelGainDiffPerCh = lStat.RunMeanRelGainDiffPerCh / lStat.CounterRelGainDiffPerCh;
      if(lStat.RunRmsRelGainDiffPerCh / lStat.CounterRelGainDiffPerCh - lStat.RunMeanRelGainDiffPerCh*lStat.RunMeanRelGainDiffPerCh){
       lStat.RunRmsRelGainDiffPerCh = sqrt( lStat.RunRmsRelGainDiffPerCh / lStat.CounterRelGainDiffPerCh - lStat.RunMeanRelGainDiffPerCh*lStat.RunMeanRelGainDiffPerCh );
      }else{
       lStat.RunRmsRelGainDiffPerCh = 0;
      }
    }else{
      lStat.RunMeanRelGainDiffPerCh = 0;
      lStat.RunRmsRelGainDiffPerCh = 0;
    }

    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(0),lDetId, lStat.RunMeanRatioPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(1),lDetId, lStat.RunMeanGainPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(2),lDetId, lStat.RunMeanDBGainPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(3),lDetId, lStat.RunMeanGainDiffPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(4),lDetId, lStat.RunMeanRelGainDiffPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(5),lDetId, lStat.RunMeanZeroLightPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(6),lDetId, lStat.RunMeanTickHeightPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(7),lDetId, lStat.RunRmsRatioPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(8),lDetId, lStat.RunRmsGainPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(9),lDetId, lStat.RunRmsDBGainPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(10),lDetId, lStat.RunRmsGainDiffPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(11),lDetId, lStat.RunRmsRelGainDiffPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(12),lDetId, lStat.RunRmsZeroLightPerCh);
    histManager_.fillTkHistoMap( histManager_.tkHistoMapPointer(13),lDetId, lStat.RunRmsTickHeightPerCh);

  }


  LogInfo("SiStripSpyGainAnalysis") << "WriteDQMStore ? " << writeDQMStore_ 
				     << " file name = " << dqmStoreFileName_
				     << std::endl;
  if (writeDQMStore_) dqm_->save(dqmStoreFileName_);

}



//
// Define as a plug-in
//

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SiStripSpyGainAnalysis);


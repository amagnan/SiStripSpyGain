#ifndef UserCode_SiStripSpyGain_GMHistograms_HH
#define UserCode_SiStripSpyGain_GMHistograms_HH

#include <sstream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DQM/SiStripCommon/interface/TkHistoMap.h"

#include "DQM/SiStripMonitorHardware/interface/SiStripSpyUtilities.h"
#include "DQM/SiStripMonitorHardware/interface/HistogramBase.hh"

class GMHistograms: public HistogramBase {

 public:

  struct Quantities {
    double ZeroLight; // digitalLow averaged over all channels per one event
    double TickHeight; // digitalHigh averaged over all channels per one event
    double Ratio; // ratio averaged over all channels per one event
    double Gain; // gain averaged over all channels per one event
    double DBGain; // gain from DB averaged over all channels per one event
    double GainDiff; // gain from DB - gain from Spy averaged over all channels per one event
    double RelGainDiff; // (gain from DB - gain from Spy) / Gain from DB averaged over all channels per one event
  };

  struct Statistics {
    double RunMeanZeroLightPerCh; // digitalLow averaged over all events in a run per one channel
    double RunRmsZeroLightPerCh;
    double CounterZeroLightPerCh;

    double RunMeanTickHeightPerCh; // digitalHigh averaged over all events in a run per one channel
    double RunRmsTickHeightPerCh;
    double CounterTickHeightPerCh;

    double RunMeanRatioPerCh; // ratio averaged over all events in a run per one channel
    double RunRmsRatioPerCh;
    double CounterRatioPerCh;

    double RunMeanGainPerCh; // gain averaged over all events in a run per one channel
    double RunRmsGainPerCh;
    double CounterGainPerCh;

    double RunMeanDBGainPerCh; // gain from DB averaged over all events in a run per one channel
    double RunRmsDBGainPerCh;
    double CounterDBGainPerCh;

    double RunMeanGainDiffPerCh; // gain from DB - gain from Spy averaged over all events in a run per one channel
    double RunRmsGainDiffPerCh;
    double CounterGainDiffPerCh;

    double RunMeanRelGainDiffPerCh; // (gain from DB - gain from Spy) / Gain from DB averaged over all events in a run per one channel
    double RunRmsRelGainDiffPerCh;
    double CounterRelGainDiffPerCh;

  };

  GMHistograms();
  
  ~GMHistograms();
  
  //initialise histograms
  void initialise(const edm::ParameterSet& iConfig,
		  std::ostringstream* pDebugStream
		  );


  void initialiseSpecificChannel(const edm::ParameterSet& iConfig,
				 std::ostringstream* pDebugStream,
				 const unsigned aRef);

  //book the top level histograms
  void bookTopLevelHistograms(DQMStore* dqm);
  void bookChannelHistograms(unsigned aRef, unsigned aPair);

  void fillRunMeanAndRunRmsHistograms(uint16_t fedCh,uint16_t fedId, Statistics aStatisticElement);

  void fillVsTimeHistograms(const unsigned aCh, const unsigned aPair, const Quantities & aElement, const double aTime);

  void fillAllEvtAllChHistograms(const double aDBGain, const double aSpyGain);

//  std::string tkHistoMapName(unsigned int);

  bool tkHistoMapEnabled(unsigned int aIndex=0);
 
  TkHistoMap * tkHistoMapPointer(unsigned int aIndex=0);

/*
  std::string tkHistoMapName(unsigned int aIndex=0);

  TkHistoMap * tkHistoMapPointer(unsigned int aIndex=0);
*/

 protected:
  
 private:

  void fill2DHistogram(HistogramConfig histogram, 
		       double xvalue, 
		       double yvalue, 
		       double weight);

  //histos

  unsigned nChannels_;

  //vsTime
  std::vector<HistogramConfig> ZeroLightvsTime_;
  std::vector<HistogramConfig> TickHeightvsTime_;
  std::vector<HistogramConfig> RatiovsTime_;
  std::vector<HistogramConfig> GainvsTime_;
  std::vector<HistogramConfig> DBGainvsTime_;
  std::vector<HistogramConfig> GainDiffvsTime_;
  std::vector<HistogramConfig> RelGainDiffvsTime_;

  //1D
  HistogramConfig allEvtAllChGainDiff_;
  HistogramConfig allEvtAllChGainRatio_;
  HistogramConfig RunMeanGainPerCh_;
  HistogramConfig RunMeanDBGainPerCh_;
  HistogramConfig RunRmsGainPerCh_;
  HistogramConfig RunRmsDBGainPerCh_;
  HistogramConfig RunMeanGainDiffPerCh_;
  HistogramConfig RunMeanRelGainDiffPerCh_;
  HistogramConfig RunRmsGainDiffPerCh_;
  HistogramConfig RunRmsRelGainDiffPerCh_;

  //2D
  HistogramConfig allEvtAllChGainScatter_;
  HistogramConfig fedIdvsChannelIdRunMeanZeroLightPerCh_;
  HistogramConfig fedIdvsChannelIdRunMeanTickHeightPerCh_;
  HistogramConfig fedIdvsChannelIdRunMeanRatioPerCh_;
  HistogramConfig fedIdvsChannelIdRunMeanGainPerCh_;
  HistogramConfig fedIdvsChannelIdRunMeanDBGainPerCh_;
  HistogramConfig fedIdvsChannelIdRunRmsZeroLightPerCh_;
  HistogramConfig fedIdvsChannelIdRunRmsTickHeightPerCh_;
  HistogramConfig fedIdvsChannelIdRunRmsRatioPerCh_;
  HistogramConfig fedIdvsChannelIdRunRmsGainPerCh_;
  HistogramConfig fedIdvsChannelIdRunRmsDBGainPerCh_;


  std::string tkMapConfigName_;
  TkHistoMap *tkmap_[14];
  HistogramConfig tkMapConfig_;


};//class



#endif //UserCode_SiStripSpyGain_GMHistograms_HH



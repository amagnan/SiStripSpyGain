#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForHardwareSystems.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQM/SiStripMonitorHardware/interface/SiStripFEDSpyBuffer.h"
#include "UserCode/SiStripSpyGain/interface/GMHistograms.h"

using edm::LogError;
using edm::LogWarning;
using edm::LogInfo;


GMHistograms::GMHistograms()
{
  dqm_ = 0;
  nChannels_ = 0;
}

GMHistograms::~GMHistograms()
{
}
 

void 
GMHistograms::initialise(const edm::ParameterSet& iConfig,
			  std::ostringstream* pDebugStream){

  nChannels_ = (iConfig.exists("ChannelList") ? (iConfig.getUntrackedParameter<std::vector<unsigned> >("ChannelList")).size() : 0);
  ZeroLightvsTime_.resize(3*nChannels_);
  TickHeightvsTime_.resize(3*nChannels_);
  RatiovsTime_.resize(3*nChannels_);
  GainvsTime_.resize(3*nChannels_);
  DBGainvsTime_.resize(3*nChannels_);
  GainDiffvsTime_.resize(3*nChannels_);
  RelGainDiffvsTime_.resize(3*nChannels_);

  for (unsigned iCh(0); iCh<3*nChannels_; ++iCh){
    initialiseSpecificChannel(iConfig,pDebugStream,iCh);
  }

  getConfigForHistogram(allEvtAllChGainDiff_,"allEvtAllChGainDiff",iConfig,pDebugStream);
  getConfigForHistogram(allEvtAllChGainRatio_,"allEvtAllChGainRatio",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunMeanZeroLightPerCh_,"fedIdvsChannelIdRunMeanZeroLightPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunMeanTickHeightPerCh_,"fedIdvsChannelIdRunMeanTickHeightPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunMeanRatioPerCh_,"fedIdvsChannelIdRunMeanRatioPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunMeanGainPerCh_,"fedIdvsChannelIdRunMeanGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunMeanDBGainPerCh_,"fedIdvsChannelIdRunMeanDBGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunRmsZeroLightPerCh_,"fedIdvsChannelIdRunRmsZeroLightPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunRmsTickHeightPerCh_,"fedIdvsChannelIdRunRmsTickHeightPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunRmsRatioPerCh_,"fedIdvsChannelIdRunRmsRatioPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunRmsGainPerCh_,"fedIdvsChannelIdRunRmsGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(fedIdvsChannelIdRunRmsDBGainPerCh_,"fedIdvsChannelIdRunRmsDBGainPerCh",iConfig,pDebugStream);

  getConfigForHistogram(allEvtAllChGainScatter_,"allEvtAllChGainScatter",iConfig,pDebugStream);

  getConfigForHistogram(RunMeanGainPerCh_,"RunMeanGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunMeanDBGainPerCh_,"RunMeanDBGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunRmsGainPerCh_,"RunRmsGainPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunRmsDBGainPerCh_,"RunRmsDBGainPerCh",iConfig,pDebugStream);


  getConfigForHistogram(RunMeanGainDiffPerCh_,"RunMeanGainDiffPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunMeanRelGainDiffPerCh_,"RunMeanRelGainDiffPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunRmsGainDiffPerCh_,"RunRmsGainDiffPerCh",iConfig,pDebugStream);
  getConfigForHistogram(RunRmsRelGainDiffPerCh_,"RunRmsRelGainDiffPerCh",iConfig,pDebugStream);

  //tkMapConfigName_ = "TkHistoMap";
  getConfigForHistogram(tkMapConfig_,"TkHistoMap",iConfig,pDebugStream);

}

void GMHistograms::initialiseSpecificChannel(const edm::ParameterSet& iConfig,
					     std::ostringstream* pDebugStream,
					     const unsigned aRef){

  getConfigForHistogram(ZeroLightvsTime_[aRef],"ZeroLightvsTime",iConfig,pDebugStream);
  getConfigForHistogram(TickHeightvsTime_[aRef],"TickHeightvsTime",iConfig,pDebugStream);
  getConfigForHistogram(RatiovsTime_[aRef],"RatiovsTime",iConfig,pDebugStream);
  getConfigForHistogram(GainvsTime_[aRef],"GainvsTime",iConfig,pDebugStream);
  getConfigForHistogram(DBGainvsTime_[aRef],"DBGainvsTime",iConfig,pDebugStream);
  getConfigForHistogram(GainDiffvsTime_[aRef],"GainDiffvsTime",iConfig,pDebugStream);
  getConfigForHistogram(RelGainDiffvsTime_[aRef],"RelGainDiffvsTime",iConfig,pDebugStream);

}


void GMHistograms::fillRunMeanAndRunRmsHistograms(uint16_t fedCh,uint16_t fedId, Statistics aStatisticElement){
  fill2DHistogram(fedIdvsChannelIdRunMeanZeroLightPerCh_,fedCh,fedId,aStatisticElement.RunMeanZeroLightPerCh);
  fill2DHistogram(fedIdvsChannelIdRunMeanTickHeightPerCh_,fedCh,fedId,aStatisticElement.RunMeanTickHeightPerCh);
  fill2DHistogram(fedIdvsChannelIdRunMeanRatioPerCh_,fedCh,fedId,aStatisticElement.RunMeanRatioPerCh);
  fill2DHistogram(fedIdvsChannelIdRunMeanGainPerCh_,fedCh,fedId,aStatisticElement.RunMeanGainPerCh);
  fill2DHistogram(fedIdvsChannelIdRunMeanDBGainPerCh_,fedCh,fedId,aStatisticElement.RunMeanDBGainPerCh);
  fill2DHistogram(fedIdvsChannelIdRunRmsZeroLightPerCh_,fedCh,fedId,aStatisticElement.RunRmsZeroLightPerCh);
  fill2DHistogram(fedIdvsChannelIdRunRmsTickHeightPerCh_,fedCh,fedId,aStatisticElement.RunRmsTickHeightPerCh);
  fill2DHistogram(fedIdvsChannelIdRunRmsRatioPerCh_,fedCh,fedId,aStatisticElement.RunRmsRatioPerCh);
  fill2DHistogram(fedIdvsChannelIdRunRmsGainPerCh_,fedCh,fedId,aStatisticElement.RunRmsGainPerCh);
  fill2DHistogram(fedIdvsChannelIdRunRmsDBGainPerCh_,fedCh,fedId,aStatisticElement.RunRmsDBGainPerCh);
  fillHistogram(RunMeanGainPerCh_,aStatisticElement.RunMeanGainPerCh);
  fillHistogram(RunMeanDBGainPerCh_,aStatisticElement.RunMeanDBGainPerCh);
  fillHistogram(RunRmsGainPerCh_,aStatisticElement.RunRmsGainPerCh);
  fillHistogram(RunRmsDBGainPerCh_,aStatisticElement.RunRmsDBGainPerCh);
  fillHistogram(RunMeanGainDiffPerCh_,aStatisticElement.RunMeanGainDiffPerCh);
  fillHistogram(RunMeanRelGainDiffPerCh_,aStatisticElement.RunMeanRelGainDiffPerCh);
  fillHistogram(RunRmsGainDiffPerCh_,aStatisticElement.RunRmsGainDiffPerCh);
  fillHistogram(RunRmsRelGainDiffPerCh_,aStatisticElement.RunRmsRelGainDiffPerCh);

}

void GMHistograms::fillVsTimeHistograms(const unsigned aCh, const unsigned aPair, 
					const Quantities & aElement, const double aTime) {
  fillHistogram(ZeroLightvsTime_[3*aCh+aPair],aTime,aElement.ZeroLight);
  fillHistogram(TickHeightvsTime_[3*aCh+aPair],aTime,aElement.TickHeight);
  fillHistogram(RatiovsTime_[3*aCh+aPair],aTime,aElement.Ratio);
  fillHistogram(GainvsTime_[3*aCh+aPair],aTime,aElement.Gain);
  fillHistogram(DBGainvsTime_[3*aCh+aPair],aTime,aElement.DBGain);
  fillHistogram(GainDiffvsTime_[3*aCh+aPair],aTime,aElement.GainDiff);
  fillHistogram(RelGainDiffvsTime_[3*aCh+aPair],aTime,aElement.RelGainDiff);
}


void GMHistograms::fill2DHistogram(HistogramConfig histogram, 
				    double xvalue, 
				    double yvalue, 
				    double weight)
{
   if (histogram.monitorEle) histogram.monitorEle->Fill(xvalue,yvalue,weight);
}


void GMHistograms::fillAllEvtAllChHistograms(const double aDBGain, const double aSpyGain){

  fillHistogram(allEvtAllChGainDiff_,aSpyGain-aDBGain);
  fillHistogram(allEvtAllChGainRatio_,aDBGain/aSpyGain);
  fill2DHistogram(allEvtAllChGainScatter_,aDBGain,aSpyGain,1);

}
void GMHistograms::bookChannelHistograms(unsigned aRef, unsigned aPair)
{

  std::ostringstream lRef;
  lRef << "_" << aRef << "_" << aPair ;
  unsigned lIndex = 3*aRef+aPair;
  bookProfile(ZeroLightvsTime_[lIndex],("ZeroLightvsTime"+lRef.str()).c_str(),
				     ";time;ZeroLight",
				     0,
				     1024,//maximum for ZeroLight
				     "Time",
				     "ZeroLight");


  bookProfile(TickHeightvsTime_[lIndex],("TickHeightvsTime"+lRef.str()).c_str(),
				     ";time;TickHeight",
				     0,
				     1024,//maximum for TickHeight
				     "Time",
				     "TickHeight");


  bookProfile(RatiovsTime_[lIndex],("RatiovsTime"+lRef.str()).c_str(),
				     ";time;Ratio",
				     0,
				     1024,//maximum for Ratio
				     "Time",
				     "Ratio");

  bookProfile(GainvsTime_[lIndex],("GainvsTime"+lRef.str()).c_str(),
				     ";time;Gain",
				     0,
				     1024,//maximum for Gain
				     "Time",
				     "Gain");

  bookProfile(DBGainvsTime_[lIndex],("DBGainvsTime"+lRef.str()).c_str(),
				     ";time;DBGain",
				     0,
				     1024,//maximum for DBGain
				     "Time",
				     "DBGain");


  bookProfile(GainDiffvsTime_[lIndex],("GainDiffvsTime"+lRef.str()).c_str(),
				     ";time;GainDiff",
				     0,
				     1024,//maximum for GainDiff
				     "Time",
				     "GainDiff");


  bookProfile(RelGainDiffvsTime_[lIndex],("RelGainDiffvsTime"+lRef.str()).c_str(),
				     ";time;RelGainDiff",
				     0,
				     1024,//maximum for RelGainDiff
				     "Time",
				     "RelGainDiff");

}

void GMHistograms::bookTopLevelHistograms(DQMStore* dqm)
{

  dqm_ = dqm;

  /*
  for (std::map<std::string,HistogramConfig>::iterator iC = histogramConfig_.begin(); iC != histogramConfig_.end(); iC++){

    LogDebug("GMHistograms") << " -- Config name : " << iC->first << ", isEnabled = " << iC->second.enabled << std::endl;
  }
  */

  //book histos
  for (unsigned iCh(0); iCh<nChannels_; ++iCh){
    for (unsigned iP(0); iP<3; ++iP){
      bookChannelHistograms(iCh,iP);
    }
  }

  //1D histo
  bookHistogram(allEvtAllChGainDiff_,"allEvtAllChGainDiff",
				       "GainDiff",
				         200, -2, 2,
                                        "SpyGain - DBGain");

  bookHistogram(allEvtAllChGainRatio_,"allEvtAllChGainRatio",
				       "GainRatio",
				         100, 0, 2,
                                        "DBGain / SpyGain");


  bookHistogram(RunMeanGainPerCh_,"RunMeanGainPerCh",
                                        "Gain",
                                        100,0,2,
					"SpyGain");

  bookHistogram(RunMeanDBGainPerCh_,"RunMeanDBGainPerCh",
                                        "DBGain",
                                        100,0,2,
					"DBGain");

  bookHistogram(RunRmsGainPerCh_,"RunRmsGainPerCh",
                                        "GainRms",
                                        100,0,0.5,
					"SpyGainRms");

  bookHistogram(RunRmsDBGainPerCh_,"RunRmsDBGainPerCh",
                                        "DBGainRms",
                                        100,0,0.5,
					"DBGainRms");

  bookHistogram(RunMeanGainDiffPerCh_,"RunMeanGainDiffPerCh",
                                        "GainDiff",
                                        100,-2,2,
					"SpyGain - DBGain");

  bookHistogram(RunMeanRelGainDiffPerCh_,"RunMeanRelGainDiffPerCh",
                                        "RelGainDiff",
                                        100,-2,2,
					"(SpyGain - DBGain)/ DBGain");

  bookHistogram(RunRmsGainDiffPerCh_,"RunRmsGainDiffPerCh",
                                        "RmsGainDiff",
                                        100,0,0.5,
					"Rms(SpyGain - DBGain)");

  bookHistogram(RunRmsRelGainDiffPerCh_,"RunRmsRelGainDiffPerCh",
                                        "RmsRelGainDiff",
                                        100,0,0.5,
					"Rms((SpyGain - DBGain)/ DBGain)");


  //2D histo
  book2DHistogram(allEvtAllChGainScatter_,"allEvtAllChGainScatter",
				       "GainScatter",
					 1000, 0, 2, 1000, 0, 2,
					 "DBGain","SpyGain");

  book2DHistogram(fedIdvsChannelIdRunMeanZeroLightPerCh_,
				      "fedIdvsChannelIdRunMeanZeroLightPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");


  book2DHistogram(fedIdvsChannelIdRunMeanTickHeightPerCh_,
				      "fedIdvsChannelIdRunMeanTickHeightPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunMeanRatioPerCh_,
				      "fedIdvsChannelIdRunMeanRatioPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunMeanGainPerCh_,
				      "fedIdvsChannelIdRunMeanGainPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunMeanDBGainPerCh_,
				      "fedIdvsChannelIdRunMeanDBGainPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunRmsZeroLightPerCh_,
				      "fedIdvsChannelIdRunRmsZeroLightPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunRmsTickHeightPerCh_,
				      "fedIdvsChannelIdRunRmsTickHeightPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunRmsRatioPerCh_,
				      "fedIdvsChannelIdRunRmsRatioPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunRmsGainPerCh_,
				      "fedIdvsChannelIdRunRmsGainPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  book2DHistogram(fedIdvsChannelIdRunRmsDBGainPerCh_,
				      "fedIdvsChannelIdRunRmsDBGainPerCh",
				      ";chId;fedId",
				      96,0,96,440,50,490,
				      "chId","fedId");

  LogInfo("GMHistograms") << " ---- Toplevel histograms book ! Number of MEs : " 
			  << dqm_->getMEs().size()
			  << std::endl;
  
  //book map after, as it creates a new folder...
  //  if (histogramConfig_[tkMapConfigName_].enabled){
  if (tkMapConfig_.enabled){
    //const std::string dqmPath = dqm_->pwd();
    tkmap_[0] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanRatioPerCh",0.,1024);
    tkmap_[1] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanGainPerCh",0.,1024);
    tkmap_[2] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanDBGainPerCh",0.,1024);
    tkmap_[3] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanGainDiffPerCh",0.,1024);
    tkmap_[4] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanRelGainDiffPerCh",0.,1024);
    tkmap_[5] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanZeroLightPerCh",0.,1024);
    tkmap_[6] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunMeanTickHeightPerCh",0.,1024);
    tkmap_[7] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsRatioPerCh",0.,1024);
    tkmap_[8] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsGainPerCh",0.,1024);
    tkmap_[9] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsDBGainPerCh",0.,1024);
    tkmap_[10] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsGainDiffPerCh",0.,1024);
    tkmap_[11] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsRelGainDiffPerCh",0.,1024);
    tkmap_[12] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsZeroLightPerCh",0.,1024);
    tkmap_[13] = new TkHistoMap("SiStrip/TkHisto","TkHMap_RunRmsTickHeightPerCh",0.,1024);
  }
  else {
    for (unsigned int iT(0); iT<14; iT++){
      tkmap_[iT] = 0;
    }
  }


}

/*
std::string GMHistograms::tkHistoMapName(unsigned int aIndex){
  return tkMapConfigName_;
}
*/
bool GMHistograms::tkHistoMapEnabled(unsigned int aIndex){
  return tkMapConfig_.enabled;
}

TkHistoMap * GMHistograms::tkHistoMapPointer(unsigned int aIndex){
  assert(aIndex < 14);
  return tkmap_[aIndex];
}


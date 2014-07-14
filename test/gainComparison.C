#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLine.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <iomanip>

int gainComparison() {

  const unsigned nRuns = 2;

  TFile *fLoic[nRuns];
  fLoic[0] = TFile::Open("Gains_195915_to_199021.root");
  fLoic[1] = TFile::Open("Gains_202093_to_203742.root");


  TFile *fin[nRuns];
  fin[0] = TFile::Open("DQMStore_merged_new.root");
  fin[1] = TFile::Open("DQMStore_202209.root");

  TString runnum[nRuns] = {"198374","202209"};

  const unsigned nIds = 51;
  const unsigned nDets = nIds*3;
  double spyGain[nRuns][nDets];
  double spyGainErr[nRuns][nDets];
  double spyHigh[nRuns][nDets];
  double spyHighErr[nRuns][nDets];
  double spyLow[nRuns][nDets];
  double spyLowErr[nRuns][nDets];
  bool valid[nRuns][nDets];
  double calibGain[nRuns][nDets*2];
  //double calibGainErr[nRuns][nDets*2];


  unsigned detids[nIds] = {
    470079220,369153672,369153636,369153668,369153640,
    369153676,369153660,369153644,369153656,402666258,
    402666261,402666257,402666262,369173524,369173540,
    369173556,369173560,369173544,369173528,369173532,
    369173548,369173564,369124374,369124381,369124373,
    369124378,369124377,369124382,436316292,436316296,
    436316300,436316304,436316308,436316312,436316316,
    436316320,436316324,436316328,436316332,436232454,
    436232466,436232469,436232474,436232457,436232462,
    436232465,436232470,436232473,436232453,436232458,
    436232461};

  //std::vector<unsigned> detidVec;

  for (unsigned iR(0); iR<nRuns; ++iR){//runs

    std::cout << " - Processing run " << runnum[iR] << std::endl;
    fin[iR]->cd("DQMData/SiStrip/ReadoutView/SpyMonitoringSummary/");

    for (unsigned int iD(0); iD<nIds; ++iD){//detids
      
      std::cout << " - Processing detid " << detids[iD] << std::endl;

      //detidVec.push_back(detids[iD]);

      for (unsigned iP(0); iP<3; ++iP){//pairs
	unsigned i = 3*iD+iP;
	std::ostringstream lName;
	TProfile *lProf = 0;
	valid[iR][i] = false;

	lName.str("");
	lName << "GainvsTime_" << iD << "_" << iP;
	lProf  = (TProfile*)gDirectory->Get(lName.str().c_str());
	if (lProf->GetEntries()==0) continue;
	spyGain[iR][i] = lProf->GetMean(2);
	spyGainErr[iR][i] = lProf->GetMeanError(2);
	valid[iR][i] = true;
	
	lName.str("");
	lName << "TickHeightvsTime_" << iD << "_" << iP;
	lProf  = (TProfile*)gDirectory->Get(lName.str().c_str());
	if (lProf->GetEntries()==0) continue;
	spyHigh[iR][i] = lProf->GetMean(2);
	spyHighErr[iR][i] = lProf->GetMeanError(2);
	
	lName.str("");
	lName << "ZeroLightvsTime_" << iD << "_" << iP;
	lProf  = (TProfile*)gDirectory->Get(lName.str().c_str());
	if (lProf->GetEntries()==0) continue;
	spyLow[iR][i] = lProf->GetMean(2);
	spyLowErr[iR][i] = lProf->GetMeanError(2);
	
      }//pairs
      
    }//detids

    fLoic[iR]->cd("SiStripCalib");

    TTree *t = (TTree*)gDirectory->Get("APVGain");
    unsigned detid = 0;
    double gain = 0;
    unsigned apvidx = 0;

    t->SetBranchAddress("DetId",&detid);
    t->SetBranchAddress("APVId",&apvidx);
    t->SetBranchAddress("Gain",&gain);

    std::cout << "--- Processing tree with: " << t->GetEntries() << " entries." << std::endl;
    for (Long64_t ievt=0; ievt<t->GetEntries();ievt++) {//apvs
      
      if (ievt%1000 == 0) std::cout << "--- ... Processing entry: " << ievt << std::endl;
      
      t->GetEntry(ievt);

      for (unsigned iD(0); iD<nIds; ++iD){//detids
	if (detid == detids[iD]) {
	  unsigned i = 6*iD + apvidx;
	  calibGain[iR][i] = gain;
	  break;
	}
      }//detids

    }//apvs

  }//runs

  TCanvas *myc = new TCanvas("myc","",1000,600);
  myc->Divide(4,2);
  const unsigned nLayers = 7;
  unsigned index[nLayers+1] = {0,1,9,13,22,28,39,nIds};
  TString layer[nLayers] = {"TEC.2.6.3.2.3.75","TID.2.4.2.1.1","TIB.2.3.6.1","TIB.1.4.1.1.1.*","TIB.1.1.1.1.1.*","TOB.1.6.1.1.1.*","TOB.1.1.1.1.1.*"};

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  const unsigned nGains = 5;
  TH1F *diffVsLayer[nLayers][nGains];
 
  std::string suffix[nGains] = {"calibAPV0","calibAPV1","spyGain","spyHigh","spyLow"};

  std::cout << " ----------------------------------------------------------------------------" << std::endl
	    << " | Detid | APV pair | Run | Calib gain | Spy gain | TickHeight | ZeroLight |" << std::endl;


  TLatex lat;
  double maxVal[nLayers];
  double minVal[nLayers];

  for (unsigned iL(0); iL<nLayers; ++iL){
  maxVal[iL] = -100;
  minVal[iL] = 100;
    for (unsigned iG(0); iG<nGains; ++iG){
      std::ostringstream histName;
      histName << "diffVsLayer_" << iL << "_" << suffix[iG];
      const unsigned nBins = (index[iL+1]-index[iL])*3;
      diffVsLayer[iL][iG] = new TH1F(histName.str().c_str(),";detid;diff(%)",nBins,0,nBins);
    }

    unsigned localIndex = 0;
    for (unsigned iD(index[iL]); iD<index[iL+1]; ++iD){//detids
      for (unsigned iP(0); iP<3; ++iP){//pairs
	unsigned iC = 6*iD+2*iP;
	unsigned iS = 3*iD+iP;
	unsigned iLocal = 3*localIndex + iP;
	for (unsigned iR(0); iR<nRuns;++iR) {//runs
	  if (valid[iR][iS] == false) continue;
	  std::cout << "| " << detids[iD] << " | "
		    << iP << " | "
		    << runnum[iR] << " | " 
		    << std::setprecision(4) 
		    << calibGain[iR][iC] << ","
		    << calibGain[iR][iC+1] << " | "
		    << spyGain[iR][iS] 
	    //  << " +/- " << spyGainErr[iR][iS] 
		    << " | "
		    << spyHigh[iR][iS] 
	    //  << " +/- " << spyHighErr[iR][iS] 
		    << " | "
		    << spyLow[iR][iS]  
	    //<< " +/- " << spyLowErr[iR][iS] 
		    << "| \n" ;//<< std::endl;
	}//runs
	if (valid[0][iS] == false || valid[1][iS] == false) continue;

	double diff[5] = {
	  (calibGain[0][iC]-calibGain[1][iC])/calibGain[1][iC]*100,
	  (calibGain[0][iC+1]-calibGain[1][iC+1])/calibGain[1][iC+1]*100,
	  (spyGain[0][iS]-spyGain[1][iS])/spyGain[1][iS]*100,
	  (spyHigh[0][iS]-spyHigh[1][iS])/spyHigh[1][iS]*100,
	  (spyLow[0][iS]-spyLow[1][iS])/spyLow[1][iS]*100
	};

	std::cout << "| " << detids[iD] << " | "
		  << iP << " | "
		  << "%diff" << " | " 
		  << std::setprecision(3) 
		  << diff[0] << ","
		  << diff[1] << " | "
		  << diff[2]
	  //  << " +/- " << spyGainErr[iR][iS] 
		  << " | "
		  << diff[3]
	  //  << " +/- " << spyHighErr[iR][iS] 
		  << " | "
		  << diff[4]
	  //<< " +/- " << spyLowErr[iR][iS] 
		  << "| \n" ;//<< std::endl;
	
	//if (iD >= index[iL] && iD < index[iL+1] ) {
	diffVsLayer[iL][0]->Fill(iLocal,diff[0]);
	diffVsLayer[iL][1]->Fill(iLocal,diff[1]);
	diffVsLayer[iL][2]->Fill(iLocal,diff[2]);
	diffVsLayer[iL][3]->Fill(iLocal,diff[3]);
	diffVsLayer[iL][4]->Fill(iLocal,diff[4]);
	//}
	for (unsigned iG(0); iG<nGains; ++iG){//gains
	  if (diff[iG] > maxVal[iL]) maxVal[iL] = diff[iG];
	  if (diff[iG] < minVal[iL]) minVal[iL] = diff[iG];
	}

      }//pairs
      localIndex++; 
    }//detids
  }//layers

  unsigned drawIndex[nLayers] = {7,6,5,2,1,4,3};
  TLegend *leg = new TLegend(0.01,0.01,0.99,0.99);
  leg->SetFillColor(10);
  for (unsigned iL(0); iL<nLayers; ++iL){
    myc->cd(drawIndex[iL]);
    for (unsigned iG(0); iG<nGains; ++iG){//gains
      diffVsLayer[iL][iG]->GetYaxis()->SetRangeUser(minVal[iL],maxVal[iL]);
      diffVsLayer[iL][iG]->SetMarkerStyle(20+iG);
      diffVsLayer[iL][iG]->SetMarkerColor(iG+1);
      diffVsLayer[iL][iG]->SetLineColor(iG+1);
      if (iG==4){
	diffVsLayer[iL][iG]->SetMarkerColor(nGains+1);
	diffVsLayer[iL][iG]->SetLineColor(nGains+1);
      }

      if (iG==0) diffVsLayer[iL][iG]->Draw("PE");
      else diffVsLayer[iL][iG]->Draw("PEsame");

      if (iL==0) leg->AddEntry(diffVsLayer[iL][iG],suffix[iG].c_str(),"P");

      std::cout << layer[iL] << " " 
		<< suffix[iG] << " " 
		<< diffVsLayer[iL][iG]->GetNbinsX() << " " 
		<< diffVsLayer[iL][iG]->GetEntries() 
		<< std::endl;

    }//gains
    lat.DrawLatex(1,maxVal[iL]+0.1,layer[iL]);

    

  }//layers

  myc->cd(8);
  leg->Draw();

  myc->Update();
  myc->Print("CalibvsSpy.pdf");

  return 0;

}//main

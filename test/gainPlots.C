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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int gainPlots() {


  //TFile *fin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/comm_tracker/Strip/SpyChannel/198372/DQMStore.root");
  //TFile *fin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/comm_tracker/Strip/SpyChannel/198374/DQMStore.root");
  //TFile *fin = TFile::Open("DQMStore_202209.root");
  TFile *fin = TFile::Open("DQMStore_merged_new.root");
  //TFile *fin = TFile::Open("DQMStore_GR198374.root");

  fin->cd("DQMData/SiStrip/ReadoutView/SpyMonitoringSummary/");

  TCanvas *myc = new TCanvas("myc","",1);

  bool doZoom = true;
  bool addDBGain = false;

  gStyle->SetOptStat(0);

  const unsigned nIds = 51;
  const unsigned nDets = nIds*3;

  TProfile *lProf[nDets];

  double DBGain[nDets];
  double spyGain[nDets];
  double DBGainErr[nDets];
  double spyGainErr[nDets];

  double xDet[nDets];
  double xDetErr[nDets];


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

  for (unsigned int iD(1); iD<nIds; ++iD){//detids

    std::cout << " - Processing detid " << detids[iD] << std::endl;
    TLegend *legend = new TLegend(0.7,0.8,1,1);
    legend->SetFillColor(10);

    double lMin = 2;
    double lMax = 0;

    myc->Clear();
    myc->cd();

    //get histos min and max
    for (unsigned iP(0); iP<3; ++iP){//pairs
      std::ostringstream lName;
      lName.str("");
      lName << "GainvsTime_" << iD << "_" << iP;
      unsigned i = 3*iD+iP;
      lProf[i]  = (TProfile*)gDirectory->Get(lName.str().c_str());

      lName.str("");
      lName << "DBGainvsTime_" << iD << "_" << iP;
      TProfile *lTmp = (TProfile*)gDirectory->Get(lName.str().c_str());
      DBGain[i] = lTmp->GetMean(2);
      DBGainErr[i] = 0;

      if (lProf[i]->GetEntries()==0) continue;
      lProf[i]->SetLineColor(iP+1);
      lProf[i]->SetMarkerColor(iP+1);
      lProf[i]->SetMarkerStyle(20+iP);

      std::cout << " -- Pair " << iP << " DB=" << DBGain[i] << " <y> = " << lProf[i]->GetMean(2) << " \\pm " << lProf[i]->GetMeanError(2) << std::endl;
      
      spyGain[i] = lProf[i]->GetMean(2);
      spyGainErr[i] = lProf[i]->GetMeanError(2);

      xDet[i] = i;
      xDetErr[i]=0;

      if (lProf[i]->GetMean(2) > lMax){
	lMax = lProf[i]->GetMean(2);
      }
      if (lProf[i]->GetMean(2) < lMin){
	lMin = lProf[i]->GetMean(2);
      }
      if (addDBGain){
	if (DBGain[i]>lMax){
	  lMax = DBGain[i];
	}
	if (DBGain[i]<lMin){
	  lMin = DBGain[i];
	}
      }
    }//pairs
    //draw histos
    TF1 *func[3]; 

    for (unsigned iP(0); iP<3; ++iP){//pairs
      std::ostringstream funcName;
      funcName << "func" << iP;
      func[iP] = new TF1(funcName.str().c_str(),"[0]",0,5000);
      unsigned i = 3*iD+iP;
      if (lProf[i]->GetEntries()==0) continue;
      //lProf[i]->GetXaxis()->SetRangeUser(1545,4000);
      lProf[i]->GetYaxis()->SetRangeUser(lMin-0.1,lMax+0.1);

      if (doZoom) lProf[i]->GetXaxis()->SetRangeUser(1300,2000);

      if (iP==0) lProf[i]->Draw("P");
      else lProf[i]->Draw("Psame");
      
      if (addDBGain){
	func[iP]->SetParameter(0,DBGain[i]);
	func[iP]->SetLineColor(iP+1);
	func[iP]->Draw("same");
      }

      std::ostringstream leg;
      leg << "APV pair " << iP;
      legend->AddEntry(lProf[i],leg.str().c_str(),"PL");

    }//pairs

    legend->Draw("same");

    TLatex lat;
    std::ostringstream latex;
    latex << "Detid " << detids[iD];
    if (!doZoom) lat.DrawLatex(100,lMin-0.09,latex.str().c_str());
    else lat.DrawLatex(1350,lMin-0.09,latex.str().c_str());

    myc->Update();
    std::ostringstream lPrint;
    lPrint << "Gain_" << detids[iD];
    if (doZoom) lPrint << "_zoom";
    myc->Print(("PLOTS_198374/pdf/"+lPrint.str()+".pdf").c_str());
    myc->Print(("PLOTS_198374/png/"+lPrint.str()+".png").c_str());

    //if (iD==14) return 0;
    
  }//detids

  TH2F *gainVsDB = (TH2F*)gDirectory->Get("allEvtAllChGainScatter");
  
  gainVsDB->RebinX(10);
  gainVsDB->RebinY(10);
  gainVsDB->SetTitle("");
  gainVsDB->GetXaxis()->SetRangeUser(0.5,1.5);
  gainVsDB->GetYaxis()->SetRangeUser(0.5,1.5);
  
  //gStyle->SetOptFit(1111);
  myc->SetGridx(1);
  myc->SetGridy(1);

  gainVsDB->Draw("colz");
  TF1 *fx = new TF1("fx","x",0.5,1.5);
  fx->SetLineColor(1);
  fx->Draw("same");
  TF1 *fxp5 = new TF1("fxp5","x+0.05",0.5,1.5);
  fxp5->SetLineColor(1);
  //fxp5->Draw("same");
  TF1 *fxm5 = new TF1("fxm5","x-0.05",0.5,1.5);
  fxm5->SetLineColor(1);
  //fxm5->Draw("same");

  TF1 *pol1 = new TF1("pol1","pol1",0,2);

  gainVsDB->Fit("pol1","0");
  
  pol1->SetLineColor(6);
  pol1->Draw("same");


  TLatex lat;
  lat.SetTextColor(6);
  std::ostringstream formulaStr;
  formulaStr << "spyGain = a + b #times DBGain"; 
  lat.DrawLatex(0.6,1.54,formulaStr.str().c_str());
  formulaStr.str("");
  formulaStr << "a=" << static_cast<int>(pol1->GetParameter(0)*100+0.5)/100. << " #pm " << static_cast<int>(pol1->GetParError(0)*100+0.5)/100.;
  lat.DrawLatex(0.6,1.45,formulaStr.str().c_str());
  formulaStr.str("");
  formulaStr << "b=" << static_cast<int>(pol1->GetParameter(1)*100+0.5)/100. << " #pm " << static_cast<int>(pol1->GetParError(1)*100+0.5)/100.;
  lat.DrawLatex(0.6,1.35,formulaStr.str().c_str());

  myc->Update();
  myc->Print("PLOTS_198374/pdf/SpyvsDBGain.pdf");
  myc->Print("PLOTS_198374/png/SpyvsDBGain.png");

  //pol1->SetParameters(690./640*pol1->GetParameter(0),690./640*pol1->GetParameter(1));
  //pol1->SetLineColor(3);
  //pol1->Draw("same");

  myc->SetGridx(0);
  myc->SetGridy(0);

  TGraphErrors *grSpy = new TGraphErrors(nDets,xDet,spyGain,xDetErr,spyGainErr);
  TGraphErrors *grDB = new TGraphErrors(nDets,xDet,DBGain,xDetErr,DBGainErr);
  
  grSpy->GetXaxis()->SetTitle("APV pair arb. index");
  grSpy->SetTitle("");
  grSpy->GetYaxis()->SetTitle("Gain value");
  grSpy->SetMarkerStyle(29);
  grSpy->SetMarkerColor(1);
  grSpy->SetFillColor(1);
  grSpy->SetLineColor(1);
  grSpy->SetMinimum(0.7);
  grSpy->SetMaximum(1.3);
  grSpy->Draw("AP");

  grDB->SetFillColor(7);
  //grDB->SetLineColor(2);
  grDB->SetFillStyle(1001);

  grDB->Draw("Bsame");
  grSpy->Draw("Psame");

  TLegend *leg = new TLegend(0.7,0.8,0.9,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(grSpy,"Spy gain","P");
  leg->AddEntry(grDB,"DB gain","F");
  leg->Draw("same");

  myc->Print("PLOTS_198374/pdf/TIBTIDSpyDBGain.pdf");
  myc->Print("PLOTS_198374/png/TIBTIDSpyDBGain.png");



  return 0;
}

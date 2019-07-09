#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void plotMassPeaks(std::string Zee, std::string Zmumu){
  //Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS_EE[nBins]; 
  TH1D * massPeakSS_EE[nBins]; 
  TH1D * massPeakOS_MuMu[nBins]; 
  TH1D * massPeakSS_MuMu[nBins]; 
  TH1D * massPeakOS_MuMu_withEff[nBins]; 
  TH1D * massPeakSS_MuMu_withEff[nBins]; 
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZmumuFile = TFile::Open(Zmumu.c_str(),"read");
  
  for(int i = 0; i<nBins; i++){
    massPeakOS_EE[i] = (TH1D*)ZeeFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_EE[i] = (TH1D*)ZeeFile->Get(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_MuMu[i] = (TH1D*)ZmumuFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_MuMu[i] = (TH1D*)ZmumuFile->Get(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_MuMu_withEff[i] = (TH1D*)ZmumuFile->Get(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_MuMu_withEff[i] = (TH1D*)ZmumuFile->Get(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    //subtract the 2
    massPeakOS_MuMu_withEff[i]->Add(massPeakSS_MuMu_withEff[i],-1);


    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetBottomMargin(0.2);
    massPeakOS_MuMu[i]->GetYaxis()->SetRangeUser(0, massPeakOS_MuMu[i]->GetMaximum()*1.4);  
    massPeakOS_MuMu[i]->GetXaxis()->SetTitle("m_{ll} (GeV)");
    massPeakOS_MuMu[i]->GetYaxis()->SetTitle("Counts");
    massPeakOS_MuMu[i]->GetXaxis()->CenterTitle();
    massPeakOS_MuMu[i]->GetYaxis()->CenterTitle();

    massPeakOS_MuMu[i]->SetStats(0);
    
    massPeakOS_MuMu[i]->SetLineColor(kBlack);
    massPeakOS_MuMu[i]->SetFillColor(kOrange-2);
    massPeakOS_MuMu[i]->Draw("HIST");
    massPeakOS_EE[i]->SetLineColor(kBlack);
    massPeakOS_EE[i]->SetFillColor(kBlack);
    massPeakOS_EE[i]->SetFillStyle(3345);
    massPeakOS_EE[i]->Draw("HIST same");
    
    massPeakSS_MuMu[i]->SetMarkerStyle(8);
    massPeakSS_MuMu[i]->SetMarkerColor(kRed);

    massPeakSS_EE[i]->SetMarkerStyle(21);
    massPeakSS_EE[i]->SetMarkerColor(kBlue);
    
    massPeakSS_EE[i]->Draw("p same");
    massPeakSS_MuMu[i]->Draw("p same");

    TLegend * leg = new TLegend(0.28,0.72,0.63,0.87);
    leg->SetBorderSize(0);
    leg->AddEntry(massPeakOS_MuMu[i],"Z #rightarrow #mu^{+}#mu^{-}","f");
    leg->AddEntry(massPeakOS_EE[i],"Z #rightarrow e^{+}e^{-}","f");
    leg->AddEntry(massPeakSS_MuMu[i],"Z #rightarrow #mu#mu (same sign)","p");
    leg->AddEntry(massPeakSS_EE[i],"Z #rightarrow ee (same sign)","p");
    leg->Draw("same");
    
    TLegend * leg2 = new TLegend(0.65,0.72+0.15*0.25,0.8,0.87);
    leg2->SetBorderSize(0);
    leg2->AddEntry((TObject*)0,"2018 PbPb","");
    leg2->AddEntry((TObject*)0,Form("%d-%d %%",c.getCentBinLow(i), c.getCentBinHigh(i)),"");
    leg2->AddEntry((TObject*)0,"p_{T}^{l} > 20 GeV","");
    leg2->Draw("same");


    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.png",c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.pdf",c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.C",c.getCentBinLow(i),c.getCentBinHigh(i)));
  
    delete leg2; 
    delete leg;
    delete c1;
  }


  
  const char * labels[7] = {"0-5%","5-10%","10-30%","30-50%", "50-70%", "70-90%", "0-100%"};
  float TAA[7] = {26.0, 20.5, 11.5, 3.82, 0.934, 0.152, 5.61};
  float Nmb = 7700 * 1606.05 * 1000.0;//glauber xsection is 7700 mb, second number is lumi, third converts from ub to mb
  float scaleFactor[7];
  for(int i = 0; i<7; i++){
    scaleFactor[i] = 1.0/TAA[i]/Nmb;
    if(i ==0 || i==1) scaleFactor[i] = scaleFactor[i] * 20.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor[i] = scaleFactor[i] * 5.0;
  }

  gStyle->SetErrorX(0);

  TH1D * yieldPlot = new TH1D("yieldPlot","",7,0,7);

  yieldPlot->SetBinContent(1,massPeakOS_MuMu_withEff[0]->Integral()*scaleFactor[0]);
  yieldPlot->SetBinError(1,TMath::Sqrt(massPeakOS_MuMu_withEff[0]->Integral())*scaleFactor[0]);
  yieldPlot->SetBinContent(2,massPeakOS_MuMu_withEff[1]->Integral()*scaleFactor[1]);
  yieldPlot->SetBinError(2,TMath::Sqrt(massPeakOS_MuMu_withEff[1]->Integral())*scaleFactor[1]);
  yieldPlot->SetBinContent(3,massPeakOS_MuMu_withEff[13]->Integral()*scaleFactor[2]);
  yieldPlot->SetBinError(3,TMath::Sqrt(massPeakOS_MuMu_withEff[13]->Integral())*scaleFactor[2]);
  yieldPlot->SetBinContent(4,massPeakOS_MuMu_withEff[14]->Integral()*scaleFactor[3]);
  yieldPlot->SetBinError(4,TMath::Sqrt(massPeakOS_MuMu_withEff[14]->Integral())*scaleFactor[3]);
  yieldPlot->SetBinContent(5,massPeakOS_MuMu_withEff[15]->Integral()*scaleFactor[4]);
  yieldPlot->SetBinError(5,TMath::Sqrt(massPeakOS_MuMu_withEff[15]->Integral())*scaleFactor[4]);
  yieldPlot->SetBinContent(6,massPeakOS_MuMu_withEff[16]->Integral()*scaleFactor[5]);
  yieldPlot->SetBinError(6,TMath::Sqrt(massPeakOS_MuMu_withEff[16]->Integral())*scaleFactor[5]);
  yieldPlot->SetBinContent(7,massPeakOS_MuMu_withEff[11]->Integral()*scaleFactor[6]);
  yieldPlot->SetBinError(7,TMath::Sqrt(massPeakOS_MuMu_withEff[11]->Integral())*scaleFactor[6]);
  for(int i = 1; i<8; i++){
    yieldPlot->GetXaxis()->SetBinLabel(i, labels[i-1]);
    yieldPlot->GetXaxis()->ChangeLabel(i,45);
  }
  
  yieldPlot->SetMarkerStyle(8);
  yieldPlot->SetMarkerColor(kBlack);
  yieldPlot->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z}");
  yieldPlot->GetYaxis()->SetRangeUser(0,yieldPlot->GetMaximum()*1.3);

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  yieldPlot->SetStats(0);
  yieldPlot->Draw();
  c1->SaveAs("plots/yields/yields.png");
  c1->SaveAs("plots/yields/yields.pdf");
  c1->SaveAs("plots/yields/yields.C");


  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: massPeakPlots <Z2EE file> <Z2mumu file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu = argv[2];
   
  plotMassPeaks(Zee, Zmumu);
  return 0; 
}

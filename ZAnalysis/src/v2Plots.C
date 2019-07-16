#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

TH1D * convertProfileToHistogram(TProfile * input, std::string title){
  int numberOfBins = input->GetXaxis()->GetNbins();
  float lowerBound = input->GetXaxis()->GetBinLowEdge(1);
  float upperBound = input->GetXaxis()->GetBinUpEdge(numberOfBins);
  TH1D * target = new TH1D(title.c_str(),title.c_str(),numberOfBins,lowerBound,upperBound);
  for(int i = 0; i<numberOfBins+2; i++){
    target->SetBinContent(i, input->GetBinContent(i));
    target->SetBinError(i, input->GetBinError(i));
  }
  return target;
}

void sqrtHist(TH1D * h){
  int numberOfBins = h->GetXaxis()->GetNbins();
  for(int i = 0; i<numberOfBins+2; i++){
    h->SetBinContent(i, TMath::Sqrt(h->GetBinContent(i)));
    h->SetBinError(i, 0.5/TMath::Sqrt(h->GetBinContent(i)) ); //   uncertainty on sqrt(x) is: 1/2 * 1/sqrt(x)
  }
}

void v2Plots(std::string Zmumu, std::string Zee){
  TH1::SetDefaultSumw2();
//  Settings s = Settings();

//  CentralityTool c = CentralityTool();
//  const int nBins = c.getNCentBins();


  //Z -> MuMu
  TH1D * v2MuMu_Num = 0; 
  TH1D * v2MuMu_Denom = 0; 
  TH1D * v2MuMu_Q1Mid = 0; 
  TH1D * v2MuMu_Q2Mid = 0; 
 
  TProfile * p_v2MuMu_Num; 
  TProfile * p_v2MuMu_Denom; 
  TProfile * p_v2MuMu_Q1Mid; 
  TProfile * p_v2MuMu_Q2Mid; 
 
  TFile * ZmumuFile = TFile::Open(Zmumu.c_str(),"read");
 
  p_v2MuMu_Num = (TProfile*)ZmumuFile->Get("v2NumVsCent");
  p_v2MuMu_Denom = (TProfile*)ZmumuFile->Get("v2DenomVsCent");
  p_v2MuMu_Q1Mid = (TProfile*)ZmumuFile->Get("v2Q1MidVsCent");
  p_v2MuMu_Q2Mid = (TProfile*)ZmumuFile->Get("v2Q2MidVsCent");

  v2MuMu_Num = convertProfileToHistogram(p_v2MuMu_Num, "h_v2NumVsCent");
  v2MuMu_Denom = convertProfileToHistogram(p_v2MuMu_Denom, "h_v2DenomVsCent");
  v2MuMu_Q1Mid = convertProfileToHistogram(p_v2MuMu_Q1Mid, "h_v2Q1MidVsCent");
  v2MuMu_Q2Mid = convertProfileToHistogram(p_v2MuMu_Q2Mid, "h_v2Q2MidVsCent");
  v2MuMu_Num->Print("All");
  
  v2MuMu_Denom->Multiply(v2MuMu_Q1Mid);
  v2MuMu_Denom->Divide(v2MuMu_Q2Mid);
  sqrtHist(v2MuMu_Denom);
  v2MuMu_Denom->Print("All");
 
  TH1D * v2 = (TH1D*)v2MuMu_Num->Clone("v2_MuMu");
  v2->Divide(v2MuMu_Denom);
  v2->Print("All");
 
  
  //Z -> EE
  TH1D * v2EE_Num = 0; 
  TH1D * v2EE_Denom = 0; 
  TH1D * v2EE_Q1Mid = 0; 
  TH1D * v2EE_Q2Mid = 0; 
 
  TProfile * p_v2EE_Num; 
  TProfile * p_v2EE_Denom; 
  TProfile * p_v2EE_Q1Mid; 
  TProfile * p_v2EE_Q2Mid; 

  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  p_v2EE_Num = (TProfile*)ZeeFile->Get("v2NumVsCent");
  p_v2EE_Denom = (TProfile*)ZeeFile->Get("v2DenomVsCent");
  p_v2EE_Q1Mid = (TProfile*)ZeeFile->Get("v2Q1MidVsCent");
  p_v2EE_Q2Mid = (TProfile*)ZeeFile->Get("v2Q2MidVsCent");

  v2EE_Num = convertProfileToHistogram(p_v2EE_Num, "h_v2NumVsCent_ee");
  v2EE_Denom = convertProfileToHistogram(p_v2EE_Denom, "h_v2DenomVsCent_ee");
  v2EE_Q1Mid = convertProfileToHistogram(p_v2EE_Q1Mid, "h_v2Q1MidVsCent_ee");
  v2EE_Q2Mid = convertProfileToHistogram(p_v2EE_Q2Mid, "h_v2Q2MidVsCent_ee");
  v2EE_Num->Print("All");
  
  v2EE_Denom->Multiply(v2EE_Q1Mid);
  v2EE_Denom->Divide(v2EE_Q2Mid);
  sqrtHist(v2EE_Denom);
  v2EE_Denom->Print("All");
 
  TH1D * v2EE = (TH1D*)v2EE_Num->Clone("v2_EE");
  v2EE->Divide(v2EE_Denom);
  v2EE->Print("All");

  TH1D * v2Combo = (TH1D*)v2->Clone("v2_Combo");
  for(int i = 0; i<v2Combo->GetXaxis()->GetNbins()+2; i++){
    float mu = v2->GetBinContent(i);
    float muErr = v2->GetBinError(i);
    float e = v2EE->GetBinContent(i);
    float eErr = v2EE->GetBinError(i);

    //calculate a weighted mean
    float point = mu/(muErr*muErr)+e/(eErr*eErr);
    float norm = 1.0/(muErr*muErr) + 1.0/(eErr*eErr);

    v2Combo->SetBinContent(i, point/norm ); 
    v2Combo->SetBinError(i, TMath::Sqrt(1.0/norm));
  }


  //***********************************************************************
  //actual creation of plot here
  TH1D * v2ComboPlot = new TH1D("v2ComboPlot","",7,0.15,7.15);
  TH1D * v2Plot = new TH1D("v2Plot","",7,0,7);
  TH1D * v2EEPlot = new TH1D("v2EEPlot","",7,-0.15,6.85);
  TH1D * ATLAS = new TH1D("ATLAS","",7,0,7);
  v2ComboPlot->SetBinContent(2,v2Combo->GetBinContent(v2Combo->FindBin(21)));
  v2ComboPlot->SetBinError(2,v2Combo->GetBinError(v2Combo->FindBin(21)));
  v2ComboPlot->SetBinContent(3,v2Combo->GetBinContent(v2Combo->FindBin(25)));
  v2ComboPlot->SetBinError(3,v2Combo->GetBinError(v2Combo->FindBin(25)));
  v2ComboPlot->SetBinContent(4,v2Combo->GetBinContent(v2Combo->FindBin(12)));
  v2ComboPlot->SetBinError(4,v2Combo->GetBinError(v2Combo->FindBin(12)));
  v2ComboPlot->SetBinContent(5,v2Combo->GetBinContent(v2Combo->FindBin(13)));
  v2ComboPlot->SetBinError(5,v2Combo->GetBinError(v2Combo->FindBin(13)));
  v2ComboPlot->SetBinContent(6,v2Combo->GetBinContent(v2Combo->FindBin(14)));
  v2ComboPlot->SetBinError(6,v2Combo->GetBinError(v2Combo->FindBin(14)));
  v2ComboPlot->SetBinContent(7,v2Combo->GetBinContent(v2Combo->FindBin(20)));
  v2ComboPlot->SetBinError(7,v2Combo->GetBinError(v2Combo->FindBin(20)));
  
  v2Plot->SetBinContent(2,v2->GetBinContent(v2->FindBin(21)));
  v2Plot->SetBinError(2,v2->GetBinError(v2->FindBin(21)));
  v2Plot->SetBinContent(3,v2->GetBinContent(v2->FindBin(25)));
  v2Plot->SetBinError(3,v2->GetBinError(v2->FindBin(25)));
  v2Plot->SetBinContent(4,v2->GetBinContent(v2->FindBin(12)));
  v2Plot->SetBinError(4,v2->GetBinError(v2->FindBin(12)));
  v2Plot->SetBinContent(5,v2->GetBinContent(v2->FindBin(13)));
  v2Plot->SetBinError(5,v2->GetBinError(v2->FindBin(13)));
  v2Plot->SetBinContent(6,v2->GetBinContent(v2->FindBin(14)));
  v2Plot->SetBinError(6,v2->GetBinError(v2->FindBin(14)));
  v2Plot->SetBinContent(7,v2->GetBinContent(v2->FindBin(20)));
  v2Plot->SetBinError(7,v2->GetBinError(v2->FindBin(20)));

  v2EEPlot->SetBinContent(2,v2EE->GetBinContent(v2EE->FindBin(21)));
  v2EEPlot->SetBinError(2,v2EE->GetBinError(v2EE->FindBin(21)));
  v2EEPlot->SetBinContent(3,v2EE->GetBinContent(v2EE->FindBin(25)));
  v2EEPlot->SetBinError(3,v2EE->GetBinError(v2EE->FindBin(25)));
  v2EEPlot->SetBinContent(4,v2EE->GetBinContent(v2EE->FindBin(12)));
  v2EEPlot->SetBinError(4,v2EE->GetBinError(v2EE->FindBin(12)));
  v2EEPlot->SetBinContent(5,v2EE->GetBinContent(v2EE->FindBin(13)));
  v2EEPlot->SetBinError(5,v2EE->GetBinError(v2EE->FindBin(13)));
  v2EEPlot->SetBinContent(6,v2EE->GetBinContent(v2EE->FindBin(14)));
  v2EEPlot->SetBinError(6,v2EE->GetBinError(v2EE->FindBin(14)));
  v2EEPlot->SetBinContent(7,v2EE->GetBinContent(v2EE->FindBin(20)));
  v2EEPlot->SetBinError(7,v2EE->GetBinError(v2EE->FindBin(20)));
   
  ATLAS->SetBinContent(1,-0.015);
  ATLAS->SetBinError(1,-0.018);

  const char * labels[7] = {"0-80%","0-80%","0-90%","0-10%", "10-30%", "30-50%", "50-100%"};

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  gStyle->SetErrorX(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  v2Plot->SetStats(0);
  v2Plot->GetYaxis()->SetRangeUser(-0.04,0.08);

  v2Plot->SetMarkerStyle(8);
  v2Plot->SetMarkerSize(1.3);
  v2Plot->SetMarkerColor(kBlue);
  v2Plot->SetLineColor(kBlue);
  v2Plot->SetLineWidth(2);
  v2Plot->GetXaxis()->CenterTitle();
  v2Plot->GetYaxis()->CenterTitle();

  v2EEPlot->SetMarkerStyle(21);
  v2EEPlot->SetMarkerColor(kRed+1);
  v2EEPlot->SetMarkerSize(1.5);
  v2EEPlot->SetLineColor(kRed+1);
  v2EEPlot->SetLineWidth(2);
  
  v2ComboPlot->SetMarkerStyle(34);
  v2ComboPlot->SetMarkerSize(1.8);
  v2ComboPlot->SetMarkerColor(kBlack);
  v2ComboPlot->SetLineColor(kBlack);
  v2ComboPlot->SetLineWidth(2);

  v2Plot->GetYaxis()->SetTitle("v_{2}");
  v2Plot->GetXaxis()->SetTitle("Centrality");

  ATLAS->SetMarkerColor(kViolet);
  ATLAS->SetLineColor(kViolet);
  ATLAS->SetMarkerStyle(21);
  ATLAS->SetLineWidth(3);
  ATLAS->SetMarkerSize(1.5);

  for(int i = 1; i<8; i++){
    v2Plot->GetXaxis()->SetBinLabel(i, labels[i-1]);
    v2Plot->GetXaxis()->ChangeLabel(i,45);
  }
 
  TLine * l = new TLine(0,0,7,0);
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);

  
  TLegend * leg = new TLegend(0.28,0.62,0.83,0.87);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb, p_{T}^{l} > 20 GeV","");
  leg->AddEntry((TObject*)0,"3 subevent SP Method","");
  leg->AddEntry(v2Plot,"#mu^{+}#mu^{-} channel","p");
  leg->AddEntry(v2EEPlot,"e^{+}e^{-} channel","p");
  leg->AddEntry(v2ComboPlot,"combined result","p");
  leg->AddEntry(ATLAS,"2012 ATLAS (#mu#mu + ee, stat error only)","p");

  v2Plot->Draw();
  l->Draw("same");
  leg->Draw("same");
  ATLAS->Draw("same");
  v2EEPlot->Draw("same");
  v2Plot->Draw("same");
  v2ComboPlot->Draw("same");
 
  c1->SaveAs("plots/v2/v2Summary.png");
  c1->SaveAs("plots/v2/v2Summary.pdf");
  c1->SaveAs("plots/v2/v2Summary.C");
 

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: v2Plots <Z2ee file> < Z2mumu file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu = argv[2];
   
  v2Plots(Zmumu, Zee);
  return 0; 
}

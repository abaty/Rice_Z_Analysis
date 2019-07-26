#include "include/combinePoints.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TBox.h"
#include "include/HistNameHelper.h"
#include "TStyle.h"
#include "TMatrixD.h"
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

void v2Plots(std::string Zmumu, std::string Zee, std::string Zmumu_syst, std::string Zee_syst){
  TH1::SetDefaultSumw2();
//  Settings s = Settings();

//  CentralityTool c = CentralityTool();
//  const int nBins = c.getNCentBins();

  CombinePoints cp = CombinePoints();

  HistNameHelper histName = HistNameHelper();

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
  TProfile * p_AvgEffMu   = (TProfile*)ZmumuFile->Get("v2AvgEffVsCent");
  p_v2MuMu_Denom = (TProfile*)ZmumuFile->Get("v2DenomVsCent");
  p_v2MuMu_Q1Mid = (TProfile*)ZmumuFile->Get("v2Q1MidVsCent");
  p_v2MuMu_Q2Mid = (TProfile*)ZmumuFile->Get("v2Q2MidVsCent");

  v2MuMu_Num = convertProfileToHistogram(p_v2MuMu_Num, "h_v2NumVsCent");
  v2MuMu_Denom = convertProfileToHistogram(p_v2MuMu_Denom, "h_v2DenomVsCent");
  v2MuMu_Q1Mid = convertProfileToHistogram(p_v2MuMu_Q1Mid, "h_v2Q1MidVsCent");
  v2MuMu_Q2Mid = convertProfileToHistogram(p_v2MuMu_Q2Mid, "h_v2Q2MidVsCent");
  v2MuMu_Num->Print("All");
  
  for(int k = 0; k<v2MuMu_Num->GetSize(); k++){
    v2MuMu_Num->SetBinContent(k, v2MuMu_Num->GetBinContent(k) / p_AvgEffMu->GetBinContent(k) );
    v2MuMu_Num->SetBinError(  k,   v2MuMu_Num->GetBinError(k) / p_AvgEffMu->GetBinContent(k) );
  }
  
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
  TProfile * p_AvgEffE   = (TProfile*)ZeeFile->Get("v2AvgEffVsCent");
  p_v2EE_Denom = (TProfile*)ZeeFile->Get("v2DenomVsCent");
  p_v2EE_Q1Mid = (TProfile*)ZeeFile->Get("v2Q1MidVsCent");
  p_v2EE_Q2Mid = (TProfile*)ZeeFile->Get("v2Q2MidVsCent");

  v2EE_Num = convertProfileToHistogram(p_v2EE_Num, "h_v2NumVsCent_ee");
  v2EE_Denom = convertProfileToHistogram(p_v2EE_Denom, "h_v2DenomVsCent_ee");
  v2EE_Q1Mid = convertProfileToHistogram(p_v2EE_Q1Mid, "h_v2Q1MidVsCent_ee");
  v2EE_Q2Mid = convertProfileToHistogram(p_v2EE_Q2Mid, "h_v2Q2MidVsCent_ee");
  v2EE_Num->Print("All");
 
  //noramlzie by average efficiency

  for(int k = 0; k<v2EE_Num->GetSize(); k++){
    v2EE_Num->SetBinContent(k, v2EE_Num->GetBinContent(k) / p_AvgEffE->GetBinContent(k) );
    v2EE_Num->SetBinError(  k,   v2EE_Num->GetBinError(k) / p_AvgEffE->GetBinContent(k) );
  }

 
  v2EE_Denom->Multiply(v2EE_Q1Mid);
  v2EE_Denom->Divide(v2EE_Q2Mid);
  sqrtHist(v2EE_Denom);
  v2EE_Denom->Print("All");
 
  TH1D * v2EE = (TH1D*)v2EE_Num->Clone("v2_EE");
  v2EE->Divide(v2EE_Denom);
  v2EE->Print("All");

//lets load our systematics

  TFile * eeSyst = TFile::Open(Zee_syst.c_str(),"read");
  TH1D * effErrorEE = (TH1D*) eeSyst->Get("effError");
  TH1D * hfErrorEE = (TH1D*) eeSyst->Get("hfError");
  TH1D * totalErrorEE = (TH1D*) eeSyst->Get("totalError");
  TFile * mumuSyst = TFile::Open(Zmumu_syst.c_str(),"read");
  TH1D * effErrorMuMu = (TH1D*) mumuSyst->Get("effError");
  TH1D * hfErrorMuMu = (TH1D*) mumuSyst->Get("hfError");
  TH1D * totalErrorMuMu = (TH1D*) mumuSyst->Get("totalError");
  
  TH1D * totalSystErrorCombo = (TH1D*) totalErrorMuMu->Clone("totalSystErrorCombo");

  TH1D * v2Combo = (TH1D*)v2->Clone("v2_Combo");
  for(int i = 0; i<v2Combo->GetXaxis()->GetNbins()+2; i++){
    float mu = v2->GetBinContent(i);
    float muStatErr = v2->GetBinError(i);
    float muEffErr = effErrorMuMu->GetBinContent(i);
    float muHfErr = hfErrorMuMu->GetBinContent(i);
    float e = v2EE->GetBinContent(i);
    float eStatErr = v2EE->GetBinError(i);
    float eEffErr = effErrorEE->GetBinContent(i);
    float eHfErr = hfErrorEE->GetBinContent(i);

    std::vector< TMatrixD > covariance;
    covariance.push_back( cp.getFullUncorrMatrix(muStatErr, eStatErr) );
    covariance.push_back( cp.getFullUncorrMatrix(muEffErr, eEffErr) );
    covariance.push_back( cp.getFullUncorrMatrix(muHfErr, eHfErr) );
    std::vector<double> combined = cp.combine(mu, e, covariance);

    //calculate a weighted mean
    v2Combo->SetBinContent(i, combined.at(0) ); 
    v2Combo->SetBinError(i, combined.at(1) );
    totalSystErrorCombo->SetBinContent(i, combined.at(2));
  }
  totalErrorMuMu->Print("All");
  totalErrorEE->Print("All");
  v2Combo->Print("All");
  totalSystErrorCombo->Print("All");
  hfErrorMuMu->Print("All");
  hfErrorEE->Print("All");


  //***********************************************************************
  //actual creation of plot here
  TH1D * v2ComboPlot = new TH1D("v2ComboPlot","",7,0.15,7.15);
  TH1D * v2Plot = new TH1D("v2Plot","",7,0,7);
  TH1D * v2EEPlot = new TH1D("v2EEPlot","",7,-0.15,6.85);
  TH1D * ATLAS = new TH1D("ATLAS","",7,0,7);

  TBox * v2ComboPlotBox[8];
  TBox * v2PlotBox[8];
  TBox * v2EEPlotBox[8];

   
  int binMapping[8] = {-1,-1,25,26,12,13,14,27}; 
  v2ComboPlot->SetBinContent(2,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[2])));
  v2ComboPlot->SetBinError(2,v2Combo->GetBinError(v2Combo->FindBin(binMapping[2])));
  v2ComboPlot->SetBinContent(3,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[3])));
  v2ComboPlot->SetBinError(3,v2Combo->GetBinError(v2Combo->FindBin(binMapping[3])));
  v2ComboPlot->SetBinContent(4,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[4])));
  v2ComboPlot->SetBinError(4,v2Combo->GetBinError(v2Combo->FindBin(binMapping[4])));
  v2ComboPlot->SetBinContent(5,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[5])));
  v2ComboPlot->SetBinError(5,v2Combo->GetBinError(v2Combo->FindBin(binMapping[5])));
  v2ComboPlot->SetBinContent(6,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[6])));
  v2ComboPlot->SetBinError(6,v2Combo->GetBinError(v2Combo->FindBin(binMapping[6])));
  v2ComboPlot->SetBinContent(7,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[7])));
  v2ComboPlot->SetBinError(7,v2Combo->GetBinError(v2Combo->FindBin(binMapping[7])));
  
  v2Plot->SetBinContent(2,v2->GetBinContent(v2->FindBin(binMapping[2])));
  v2Plot->SetBinError(2,v2->GetBinError(v2->FindBin(binMapping[2])));
  v2Plot->SetBinContent(3,v2->GetBinContent(v2->FindBin(binMapping[3])));
  v2Plot->SetBinError(3,v2->GetBinError(v2->FindBin(binMapping[3])));
  v2Plot->SetBinContent(4,v2->GetBinContent(v2->FindBin(binMapping[4])));
  v2Plot->SetBinError(4,v2->GetBinError(v2->FindBin(binMapping[4])));
  v2Plot->SetBinContent(5,v2->GetBinContent(v2->FindBin(binMapping[5])));
  v2Plot->SetBinError(5,v2->GetBinError(v2->FindBin(binMapping[5])));
  v2Plot->SetBinContent(6,v2->GetBinContent(v2->FindBin(binMapping[6])));
  v2Plot->SetBinError(6,v2->GetBinError(v2->FindBin(binMapping[6])));
  v2Plot->SetBinContent(7,v2->GetBinContent(v2->FindBin(binMapping[7])));
  v2Plot->SetBinError(7,v2->GetBinError(v2->FindBin(binMapping[7])));

  v2EEPlot->SetBinContent(2,v2EE->GetBinContent(v2EE->FindBin(binMapping[2])));
  v2EEPlot->SetBinError(2,v2EE->GetBinError(v2EE->FindBin(binMapping[2])));
  v2EEPlot->SetBinContent(3,v2EE->GetBinContent(v2EE->FindBin(binMapping[3])));
  v2EEPlot->SetBinError(3,v2EE->GetBinError(v2EE->FindBin(binMapping[3])));
  v2EEPlot->SetBinContent(4,v2EE->GetBinContent(v2EE->FindBin(binMapping[4])));
  v2EEPlot->SetBinError(4,v2EE->GetBinError(v2EE->FindBin(binMapping[4])));
  v2EEPlot->SetBinContent(5,v2EE->GetBinContent(v2EE->FindBin(binMapping[5])));
  v2EEPlot->SetBinError(5,v2EE->GetBinError(v2EE->FindBin(binMapping[5])));
  v2EEPlot->SetBinContent(6,v2EE->GetBinContent(v2EE->FindBin(binMapping[6])));
  v2EEPlot->SetBinError(6,v2EE->GetBinError(v2EE->FindBin(binMapping[6])));
  v2EEPlot->SetBinContent(7,v2EE->GetBinContent(v2EE->FindBin(binMapping[7])));
  v2EEPlot->SetBinError(7,v2EE->GetBinError(v2EE->FindBin(binMapping[7])));
   
  ATLAS->SetBinContent(1,-0.015);
  ATLAS->SetBinError(1,-0.018);
  TBox * ATLASBox;

  const char * labels[7] = {"0-80%","0-90%","20-80%","0-10%", "10-30%", "30-50%", "50-90%"};

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  gStyle->SetErrorX(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  v2Plot->SetStats(0);
  v2Plot->GetYaxis()->SetRangeUser(-0.07,0.12);

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

  
  TLegend * leg = new TLegend(0.28,0.64,0.83,0.87);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb, p_{T}^{l} > 20 GeV","");
  leg->AddEntry((TObject*)0,"3 subevent SP Method","");
  leg->AddEntry(v2Plot,"#mu^{+}#mu^{-} channel, |#eta_{#mu}|<2.4","p");
  leg->AddEntry(v2EEPlot,"e^{+}e^{-} channel, |#eta_{e}|<2.1","p");
  leg->AddEntry(v2ComboPlot,"combined result","p");
  leg->AddEntry(ATLAS,"2012 ATLAS","p");

  v2Plot->Draw();
  l->Draw("same");
  leg->Draw("same");
  ATLAS->Draw("same");
  v2EEPlot->Draw("same");
  v2Plot->Draw("same");
  v2ComboPlot->Draw("same");

  for(int i = 1; i<8; i++){
    if(i==1) histName.drawBoxAbsolute(ATLAS, i, ATLASBox, 0.014, 0.1, kViolet);
    if(i>1){
      histName.drawBoxAbsolute(v2ComboPlot, i , v2ComboPlotBox[i], totalSystErrorCombo->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i])),0.1,(Color_t)kBlack); 
      histName.drawBoxAbsolute(v2Plot, i , v2PlotBox[i],           totalErrorMuMu->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i])),0.1,(Color_t)kBlue); 
      histName.drawBoxAbsolute(v2EEPlot, i , v2EEPlotBox[i],       totalErrorEE->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i])),0.1,(Color_t)(kRed+1)); 
    }
  }
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
  if(argc != 5)
  {
    std::cout << "Usage: v2Plots <Z2ee file> < Z2mumu file> <Z2ee Syst file> <Z2mumu Syst file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu = argv[2];
  std::string Zee_syst = argv[3];
  std::string Zmumu_syst = argv[4];
   
  v2Plots(Zmumu, Zee, Zmumu_syst, Zee_syst);
  return 0; 
}

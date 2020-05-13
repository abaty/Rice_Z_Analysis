#include "TMathText.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "include/CMS_lumi.C"
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
  TH1D * v2ComboPlot = new TH1D("v2ComboPlot","",6,0.15,6.15);
  TH1D * v2Plot = new TH1D("v2Plot","",6,0,6);
  TH1D * v2EEPlot = new TH1D("v2EEPlot","",6,-0.15,5.85);
  TH1D * ATLAS = new TH1D("ATLAS","",6,0,6);

  TBox * v2ComboPlotBox[7];
  TBox * v2PlotBox[7];
  TBox * v2EEPlotBox[7];

   
  //int binMapping[8] = {-1,-1,25,12,13,14,27,26}; 
  int binMapping[7] = {-1,-1,25,12,13,14,27}; 
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
  //v2ComboPlot->SetBinContent(7,v2Combo->GetBinContent(v2Combo->FindBin(binMapping[7])));
  //v2ComboPlot->SetBinError(7,v2Combo->GetBinError(v2Combo->FindBin(binMapping[7])));
  
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
  //v2Plot->SetBinContent(7,v2->GetBinContent(v2->FindBin(binMapping[7])));
  //v2Plot->SetBinError(7,v2->GetBinError(v2->FindBin(binMapping[7])));

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
  //v2EEPlot->SetBinContent(7,v2EE->GetBinContent(v2EE->FindBin(binMapping[7])));
  //v2EEPlot->SetBinError(7,v2EE->GetBinError(v2EE->FindBin(binMapping[7])));
   
  ATLAS->SetBinContent(1,-0.015);
  ATLAS->SetBinError(1,-0.018);
  TBox * ATLASBox;

  //const char * labels[7] = {"0-80%","0-90%","0-10%", "10-30%", "30-50%", "50-90%","20-80%"};
  const char * labels[6] = {"0-80","0-90","0-10", "10-30", "30-50", "50-90"};

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  gStyle->SetErrorX(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  v2Plot->SetStats(0);
  v2Plot->GetYaxis()->SetRangeUser(-0.05,0.12);

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
  v2Plot->GetYaxis()->SetTitleOffset(1.11);
  v2Plot->GetYaxis()->SetTitleSize(0.08);
  v2Plot->GetYaxis()->SetLabelSize(0.04);
  v2Plot->GetXaxis()->SetTitle("Centrality (%)");
  v2Plot->GetXaxis()->SetTitleOffset(1.0);
  v2Plot->GetXaxis()->SetTitleSize(0.065);
  v2Plot->GetXaxis()->SetLabelSize(0.06);
  v2Plot->GetXaxis()->SetLabelOffset(0.007);

  ATLAS->SetMarkerColor(kViolet);
  ATLAS->SetLineColor(kViolet);
  ATLAS->SetMarkerStyle(21);
  ATLAS->SetLineWidth(3);
  ATLAS->SetMarkerSize(1.5);

  for(int i = 1; i<7; i++){
    v2Plot->GetXaxis()->SetBinLabel(i, labels[i-1]);
    v2Plot->GetXaxis()->ChangeLabel(i,45);
  }
 
  TLine * l = new TLine(0,0,6,0);
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);

  
  TLegend * leg = new TLegend(0.28,0.64,0.83,0.87);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb, p_{T}^{l} > 20 GeV","");
  //leg->AddEntry((TObject*)0,"3 subevent SP Method","");
  leg->AddEntry(v2Plot,"#mu^{+}#mu^{-}, |#eta_{#mu}|<2.4","p");
  leg->AddEntry(v2EEPlot,"e^{+}e^{-}, |#eta_{e}|<2.1","p");
  leg->AddEntry(v2ComboPlot,"Combined","p");
  leg->AddEntry(ATLAS,"PRL 110, 022301 (2013)","p");

  v2Plot->Draw();
  l->Draw("same");
  leg->Draw("same");
  ATLAS->Draw("same");
  v2EEPlot->Draw("same");
  v2Plot->Draw("same");
  v2ComboPlot->Draw("same");

  for(int i = 1; i<7; i++){
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
  
  c1->RedrawAxis();
  CMS_lumi(c1,0,10,1.0,true,true,false,false);
 
  c1->SaveAs("plots/v2/v2Summary.png");
  c1->SaveAs("plots/v2/v2Summary.pdf");
  c1->SaveAs("plots/v2/v2Summary.C");
 
  TCanvas * c2 = new TCanvas("c2","c2",800,500);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetBottomMargin(0.15);
  float xOff = 0.5;
  //float xOffSyst = 1;
  float centersX[5] = {5,20,40,70,95.0};
  float centersXE[5] = {2.5,17.5,37.5,67.5,92.5};
  float centersXMu[5] = {7.5,22.5,42.5,72.5,97.5};
  float xErrs[5] = {5-xOff,10-xOff,10-xOff,20-xOff,0};
  //float xSystErrs[5] = {5-xOffSyst,10-xOffSyst,10-xOffSyst,20-xOffSyst,1};
  float xSystErrs[5] = {2.,2.,2.,2.,2.};
  float centersY[5] = {0};
  float yErrs[5] = {0};
  float ySystErrs[5] = {0};
  float centersYE[5] = {0};
  float yErrsE[5] = {0};
  float ySystErrsE[5] = {0};
  float centersYMu[5] = {0};
  float yErrsMu[5] = {0};
  float ySystErrsMu[5] = {0};
  for(int i = 0; i<4; i++){
    yErrs[i] = v2ComboPlot->GetBinError(i+3);
    ySystErrs[i] = totalSystErrorCombo->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i+3]));
    centersY[i] = v2ComboPlot->GetBinContent(i+3);
    yErrsE[i] = v2EEPlot->GetBinError(i+3);
    ySystErrsE[i] = totalErrorEE->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i+3]));
    centersYE[i] = v2EEPlot->GetBinContent(i+3);
    yErrsMu[i] = v2Plot->GetBinError(i+3);
    ySystErrsMu[i] = totalErrorMuMu->GetBinContent(totalSystErrorCombo->FindBin(binMapping[i+3]));
    centersYMu[i] = v2Plot->GetBinContent(i+3);
  }
  yErrs[4] = v2ComboPlot->GetBinError(2);
  centersY[4] = v2ComboPlot->GetBinContent(2);
  yErrsE[4] = v2EEPlot->GetBinError(2);
  centersYE[4] = v2EEPlot->GetBinContent(2);
  yErrsMu[4] = v2Plot->GetBinError(2);
  centersYMu[4] = v2Plot->GetBinContent(2);

  TGraphAsymmErrors * v2ComboGraph = new TGraphAsymmErrors(5, centersX, centersY, xErrs, xErrs, yErrs, yErrs );
  TGraphAsymmErrors * v2ComboGraph_noX = new TGraphAsymmErrors(5, centersX, centersY, 0, 0, yErrs, yErrs );
  TGraphAsymmErrors * v2ComboGraphSyst = new TGraphAsymmErrors(5, centersX, centersY, xSystErrs, xSystErrs, ySystErrs, ySystErrs );
  v2ComboGraph->SetMarkerStyle(8);
  v2ComboGraph->SetMarkerSize(1.5);
  v2ComboGraph->SetMarkerColor(kBlack);
  v2ComboGraph->SetLineColor(kBlack);
  v2ComboGraph->SetLineWidth(1);  
  v2ComboGraph_noX->SetMarkerStyle(8);
  v2ComboGraph_noX->SetMarkerSize(1.);
  v2ComboGraph_noX->SetMarkerColor(kBlack);
  v2ComboGraph_noX->SetLineColor(kBlack);
  v2ComboGraph_noX->SetLineWidth(1);  
  TGraphAsymmErrors * v2EGraph = new TGraphAsymmErrors(5, centersXE, centersYE, 0, 0, yErrsE, yErrsE );
  TGraphAsymmErrors * v2EGraphSyst = new TGraphAsymmErrors(5, centersXE, centersYE, xSystErrs, xSystErrs, ySystErrsE, ySystErrsE );
  v2EGraph->SetMarkerStyle(24);
  v2EGraph->SetMarkerSize(1.);
  v2EGraph->SetMarkerColor(kRed);
  v2EGraph->SetLineColor(kRed);
  v2EGraph->SetLineWidth(1);  
  TGraphAsymmErrors * v2MuGraph = new TGraphAsymmErrors(5, centersXMu, centersYMu, 0, 0, yErrsMu, yErrsMu );
  TGraphAsymmErrors * v2MuGraphSyst = new TGraphAsymmErrors(5, centersXMu, centersYMu, xSystErrs, xSystErrs, ySystErrsMu, ySystErrsMu );
  v2MuGraph->SetMarkerStyle(25);
  v2MuGraph->SetMarkerSize(1.);
  v2MuGraph->SetMarkerColor(kBlue);
  v2MuGraph->SetLineColor(kBlue);
  v2MuGraph->SetLineWidth(1);  

  v2EGraphSyst->SetFillColor(0);
  v2MuGraphSyst->SetFillColor(0);
  v2EGraphSyst->SetLineColor(kRed);
  v2MuGraphSyst->SetLineColor(kBlue);
  v2ComboGraphSyst->SetFillStyle(0);
  v2ComboGraphSyst->SetLineColor(kBlack);
  
  TGraphAsymmErrors * ATLASg = new TGraphAsymmErrors(1);
  ATLASg->SetPoint(0,105,-0.015);
  ATLASg->SetPointError(0,0,0,0.018,0.018);
  ATLASg->SetLineColor(kViolet-1);
  ATLASg->SetMarkerColor(kViolet-1);
  ATLASg->SetMarkerStyle(21);
  ATLASg->SetMarkerSize(1.5);
  TGraphAsymmErrors * ATLASSyst = new TGraphAsymmErrors(1);
  ATLASSyst->SetPoint(0,105,-0.015);
  ATLASSyst->SetPointError(0,2.5,2.5,0.014,0.014);
  ATLASSyst->SetLineColor(kViolet-1);
  ATLASSyst->SetFillStyle(0);

  TH1D * d = new TH1D("d",";Centrality (%); v_{2}",1,0,110);
  d->SetLineColor(kWhite);
  d->GetYaxis()->SetRangeUser(-0.05,0.07);
  d->SetStats(0);
  d->GetYaxis()->CenterTitle();
  d->GetXaxis()->CenterTitle();
  d->GetYaxis()->SetTitle("v_{2}");
  d->GetYaxis()->SetTitleOffset(0.9);
  d->GetYaxis()->SetTitleSize(0.08);
  d->GetYaxis()->SetLabelSize(0.055);
  d->GetXaxis()->SetTitle("Centrality (%)");
  d->GetXaxis()->SetTitleOffset(1.0);
  d->GetXaxis()->SetTitleSize(0.065);
  d->GetXaxis()->SetLabelSize(0.06);
  d->GetXaxis()->SetLabelOffset(0.007);
  TAxis* a = d->GetXaxis();
  a->ChangeLabel(6,-1,-1,-1,-1,-1,"0-90%");
  a->SetTickLength(0);
  TGaxis * aClone = new TGaxis(0,-0.05,90,-0.05,0,90,505,"");
  aClone->SetTickLength(0.02);
  aClone->SetLabelSize(0);
  TGaxis * aClone2 = new TGaxis(0,0.07,90,0.07,0,90,505,"-");
  aClone2->SetTickLength(0.02);
  aClone2->SetLabelSize(0);

  d->Draw();

  v2ComboGraphSyst->Draw("same 5");
  ATLASSyst->Draw("same 2"); 

  TLine * l2= new TLine(0,0,110,0);
  l2->SetLineColor(kBlack);
  l2->SetLineStyle(1);
  l2->Draw("same"); 
  TLine * l3= new TLine(90,-0.05,90,0.07);
  l3->SetLineColor(kBlack);
  l3->SetLineStyle(7);
  l3->Draw("same"); 
 
  v2ComboGraph->Draw("same p");
  ATLASg->Draw("same p");

  gStyle->SetLegendBorderSize(0);
  TLegend * leg2 = new TLegend(0.3,0.57,0.7,0.89);
  leg2->SetTextSize(0.06);
  leg2->AddEntry((TObject*)0,"60 < m_{ll} < 120 GeV", "");
  leg2->AddEntry((TObject*)0,"|y_{Z}| < 2.1", "");
  leg2->AddEntry(v2ComboGraph,"Z/#gamma* #rightarrow l^{+}l^{-}", "ep");
  leg2->AddEntry(ATLASg,"PRL 110, 022301 (2013)","ep");
  leg2->Draw("same");

  aClone->Draw(); 
  aClone2->Draw(); 
 
  CMS_lumi(c2,0,10,1.5, true, true, false, false);
  c2->SaveAs("plots/v2/v2Summary_PRL.png");
  c2->SaveAs("plots/v2/v2Summary_PRL.pdf");
  c2->SaveAs("plots/v2/v2Summary_PRL.C");
  
  delete v2ComboGraph;

  d->GetYaxis()->SetRangeUser(-0.07,0.16);
  //v2ComboGraph->SetFillColor(0);
  v2EGraphSyst->SetFillStyle(0);
  v2MuGraphSyst->SetFillStyle(0);
  v2EGraphSyst->Draw("same 5");
  v2MuGraphSyst->Draw("same 5");
  v2EGraph->Draw("same p");
  v2MuGraph->Draw("same p");
  v2ComboGraph_noX->Draw("same p");

  l3->SetY1(-0.07);
  l3->SetY2(0.16);
  l3->Draw("same");
  delete leg2;
  TLegend * leg3 = new TLegend(0.4,0.57,0.8,0.87);
  leg3->SetTextSize(0.06);
  leg3->AddEntry((TObject*)0,"60 < m_{ll} < 120 GeV", "");
  leg3->AddEntry((TObject*)0,"|y_{Z}| < 2.1", "");
  leg3->AddEntry(v2ComboGraph_noX,"Z/#gamma* #rightarrow l^{+}l^{-}", "ep");
  leg3->AddEntry(ATLASg,"PRL 110, 022301 (2013)","ep");
  leg3->AddEntry(v2EGraph,"Z#rightarrow e^{+}e^{-}","ep");
  leg3->AddEntry(v2MuGraph,"Z#rightarrow #mu^{+}#mu^{-}","ep");
  leg3->SetFillStyle(0);
  leg3->Draw("same");
  CMS_lumi(c2,0,10,1.5,true,true,true,false);
  c2->SaveAs("plots/v2/v2Summary_PRL_AllChannels.png");
  c2->SaveAs("plots/v2/v2Summary_PRL_AllChannels.pdf");
  c2->SaveAs("plots/v2/v2Summary_PRL_AllChannels.C");

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

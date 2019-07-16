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
  TH1D * yields_MuMu;
  TH1D * yields_EE;
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZmumuFile = TFile::Open(Zmumu.c_str(),"read");
    
  yields_MuMu = (TH1D*)ZmumuFile->Get("yields");
  yields_EE = (TH1D*)ZeeFile->Get("yields");

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


  
  const char * labels[9] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%", "50-70%", "70-90%", "0-90%"};
  float TAA[9] = {26.0, 20.5, 14.35, 8.663, 4.976, 2.66, 0.934, 0.152, 6.2287};
  //float Nmb_mumu = 7700 * 1606.05 * 1000.0;//glauber xsection is 7700 mb, second number is lumi, third converts from ub to mb
  float Nmb_mumu = 1606.05/1618.466*9.62460e+09;//0--90%
  float scaleFactor_mumu[9];
  for(int i = 0; i<9; i++){
    scaleFactor_mumu[i] = 1.0/TAA[i]/Nmb_mumu;
    //if(i ==0 || i==1) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 20.0;
    //if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 10.0;
    //if(i == 6 || i==7) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 5.0;
    if(i ==0 || i==1) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 4.5;
  }
  //float Nmb_ee = 7700 * 1603.392 * 1000.0;//glauber xsection is 7700 mb, second number is lumi, third converts from ub to mb
  float Nmb_ee = 1603.392/1618.466*9.62460e+09;//0-90%
  float scaleFactor_ee[9];
  for(int i = 0; i<9; i++){
    scaleFactor_ee[i] = 1.0/TAA[i]/Nmb_ee;
    //if(i ==0 || i==1) scaleFactor_ee[i] = scaleFactor_ee[i] * 20.0;
    //if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_ee[i] = scaleFactor_ee[i] * 10.0;
    //if(i == 6 || i==7) scaleFactor_ee[i] = scaleFactor_ee[i] * 5.0;
    if(i ==0 || i==1) scaleFactor_ee[i] = scaleFactor_ee[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_ee[i] = scaleFactor_ee[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_ee[i] = scaleFactor_ee[i] * 4.5;
  }

  gStyle->SetErrorX(0);

  TH1D * yieldPlot_mumu = new TH1D("yieldPlot_mumu","",9,0,9);
  yieldPlot_mumu->SetBinContent(1,yields_MuMu->GetBinContent(0+1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinError(1,yields_MuMu->GetBinError(0+1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinContent(2,yields_MuMu->GetBinContent(1+1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinError(2,yields_MuMu->GetBinError(1+1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinContent(3,yields_MuMu->GetBinContent(2+1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinError(3,yields_MuMu->GetBinError(2+1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinContent(4,yields_MuMu->GetBinContent(3+1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinError(4,yields_MuMu->GetBinError(3+1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinContent(5,yields_MuMu->GetBinContent(4+1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinError(5,yields_MuMu->GetBinError(4+1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinContent(6,yields_MuMu->GetBinContent(5+1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinError(6,yields_MuMu->GetBinError(5+1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinContent(7,yields_MuMu->GetBinContent(15+1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinError(7,yields_MuMu->GetBinError(15+1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinContent(8,yields_MuMu->GetBinContent(16+1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinError(8,yields_MuMu->GetBinError(16+1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinContent(9,yields_MuMu->GetBinContent(25+1)*scaleFactor_mumu[8]);
  yieldPlot_mumu->SetBinError(9,yields_MuMu->GetBinError(25+1)*scaleFactor_mumu[8]);
  
  TH1D * yieldPlot_ee = new TH1D("yieldPlot_ee","",9,-0.25,8.75);
  yieldPlot_ee->SetBinContent(1,yields_EE->GetBinContent(0+1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinError(1,yields_EE->GetBinError(0+1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinContent(2,yields_EE->GetBinContent(1+1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinError(2,yields_EE->GetBinError(1+1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinContent(3,yields_EE->GetBinContent(2+1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinError(3,yields_EE->GetBinError(2+1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinContent(4,yields_EE->GetBinContent(3+1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinError(4,yields_EE->GetBinError(3+1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinContent(5,yields_EE->GetBinContent(4+1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinError(5,yields_EE->GetBinError(4+1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinContent(6,yields_EE->GetBinContent(5+1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinError(6,yields_EE->GetBinError(5+1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinContent(7,yields_EE->GetBinContent(15+1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinError(7,yields_EE->GetBinError(15+1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinContent(8,yields_EE->GetBinContent(16+1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinError(8,yields_EE->GetBinError(16+1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinContent(9,yields_EE->GetBinContent(25+1)*scaleFactor_ee[8]);
  yieldPlot_ee->SetBinError(9,yields_EE->GetBinError(25+1)*scaleFactor_ee[8]);

  TH1D * yieldCombo = new TH1D("yieldCombo","",9,0.25,9.25);
  for(int i = 0; i<yieldCombo->GetXaxis()->GetNbins()+2; i++){
    float mu = yieldPlot_mumu->GetBinContent(i);
    float muErr = yieldPlot_mumu->GetBinError(i);
    float e = yieldPlot_ee->GetBinContent(i);
    float eErr = yieldPlot_ee->GetBinError(i);

    //calculate a weighted mean
    float point = mu/(muErr*muErr)+e/(eErr*eErr);
    float norm = 1.0/(muErr*muErr) + 1.0/(eErr*eErr);

    yieldCombo->SetBinContent(i, point/norm ); 
    yieldCombo->SetBinError(i, TMath::Sqrt(1.0/norm));
  }


  for(int i = 1; i<10; i++){
    yieldPlot_mumu->GetXaxis()->SetBinLabel(i, labels[i-1]);
    yieldPlot_mumu->GetXaxis()->ChangeLabel(i,45);
  }
  
  yieldPlot_mumu->SetMarkerStyle(8);
  yieldPlot_mumu->SetMarkerColor(kBlue);
  yieldPlot_mumu->SetMarkerSize(1.3);
  yieldPlot_mumu->SetLineColor(kBlue);
  yieldPlot_mumu->SetLineWidth(2);
  yieldPlot_mumu->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z}");
  yieldPlot_mumu->GetYaxis()->SetRangeUser(yieldPlot_mumu->GetMinimum()*0.5,yieldPlot_mumu->GetMaximum()*1.15);

  yieldPlot_ee->SetMarkerStyle(21);
  yieldPlot_ee->SetMarkerSize(1.5);
  yieldPlot_ee->SetMarkerColor(kRed+1);
  yieldPlot_ee->SetLineColor(kRed+1);
  yieldPlot_ee->SetLineWidth(2);
  
  yieldCombo->SetMarkerStyle(34);
  yieldCombo->SetMarkerSize(1.8);
  yieldCombo->SetMarkerColor(kBlack);
  yieldCombo->SetLineColor(kBlack);
  yieldCombo->SetLineWidth(2);

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  yieldPlot_mumu->SetStats(0);
  yieldPlot_mumu->Draw();
  yieldPlot_ee->Draw("same");
  yieldCombo->Draw("same");
  
  yieldPlot_mumu->Print("All");
  yieldPlot_ee->Print("All");
  yieldCombo->Print("All");  

  TLegend * leg = new TLegend(0.21,0.22,0.81,0.43);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb, p_{T}^{l} > 20 GeV","");
  leg->AddEntry((TObject*)0,"|#eta^{l}| < 2.1","");
  leg->AddEntry(yieldPlot_mumu,"#mu^{+}#mu^{-} channel","p");
  leg->AddEntry(yieldPlot_ee,"e^{+}e^{-} channel","p");
  leg->AddEntry(yieldCombo,"combined result","p");

  leg->Draw("same");

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

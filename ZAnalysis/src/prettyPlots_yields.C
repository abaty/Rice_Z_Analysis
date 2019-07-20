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

void plotMassPeaks(std::string Zee, std::string Zmumu21, std::string Zmumu24){
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * mu24[nBins];
  TH1D * mu21[nBins];
  TH1D * e[nBins];
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZmumuFile21 = TFile::Open(Zmumu21.c_str(),"read");
  TFile * ZmumuFile24 = TFile::Open(Zmumu24.c_str(),"read");
    

  for(int i = 0; i<nBins; i++){
    mu24[i] = (TH1D*)ZmumuFile24->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    mu21[i] = (TH1D*)ZmumuFile21->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    e[i] = (TH1D*)ZeeFile->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    //scale by lumi and nMB
    mu24[i]->Scale(s.netLumi/(s.muLumi * s.Nmb));
    mu21[i]->Scale(s.netLumi/(s.muLumi * s.Nmb));
    e[i]->Scale(s.netLumi/(s.eLumi * s.Nmb));
  }


  const char * labels[10] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%", "50-70%", "70-90%", "0-90%","pp"};
  float TAA[9] = {26.0, 20.5, 14.35, 8.663, 4.976, 2.66, 0.934, 0.152, 6.2287};
  float scaleFactor_mumu[9];
  for(int i = 0; i<9; i++){
    scaleFactor_mumu[i] = 1.0/TAA[i];
    if(i ==0 || i==1) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 4.5;
  }
  float scaleFactor_ee[9];
  for(int i = 0; i<9; i++){
    scaleFactor_ee[i] = 1.0/TAA[i];
    if(i ==0 || i==1) scaleFactor_ee[i] = scaleFactor_ee[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_ee[i] = scaleFactor_ee[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_ee[i] = scaleFactor_ee[i] * 4.5;
  }

  gStyle->SetErrorX(0);

  TH1D * yieldPlot_mumu24 = new TH1D("yieldPlot_mumu24","",10,0,10);
  yieldPlot_mumu24->SetBinContent(1,mu24[0]->GetBinContent(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu24->SetBinError(1,mu24[0]->GetBinError(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu24->SetBinContent(2,mu24[1]->GetBinContent(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu24->SetBinError(2,mu24[1]->GetBinError(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu24->SetBinContent(3,mu24[2]->GetBinContent(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu24->SetBinError(3,mu24[2]->GetBinError(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu24->SetBinContent(4,mu24[3]->GetBinContent(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu24->SetBinError(4,mu24[3]->GetBinError(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu24->SetBinContent(5,mu24[4]->GetBinContent(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu24->SetBinError(5,mu24[4]->GetBinError(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu24->SetBinContent(6,mu24[5]->GetBinContent(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu24->SetBinError(6,mu24[5]->GetBinError(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu24->SetBinContent(7,mu24[15]->GetBinContent(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu24->SetBinError(7,mu24[15]->GetBinError(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu24->SetBinContent(8,mu24[16]->GetBinContent(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu24->SetBinError(8,mu24[16]->GetBinError(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu24->SetBinContent(9,mu24[25]->GetBinContent(1)*scaleFactor_mumu[8]);
  yieldPlot_mumu24->SetBinError(9,mu24[25]->GetBinError(1)*scaleFactor_mumu[8]);
  
  TH1D * yieldPlot_mumu = new TH1D("yieldPlot_mumu","",10,0,10);
  yieldPlot_mumu->SetBinContent(1,mu21[0]->GetBinContent(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinError(1,mu21[0]->GetBinError(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinContent(2,mu21[1]->GetBinContent(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinError(2,mu21[1]->GetBinError(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinContent(3,mu21[2]->GetBinContent(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinError(3,mu21[2]->GetBinError(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinContent(4,mu21[3]->GetBinContent(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinError(4,mu21[3]->GetBinError(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinContent(5,mu21[4]->GetBinContent(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinError(5,mu21[4]->GetBinError(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinContent(6,mu21[5]->GetBinContent(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinError(6,mu21[5]->GetBinError(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinContent(7,mu21[15]->GetBinContent(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinError(7,mu21[15]->GetBinError(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinContent(8,mu21[16]->GetBinContent(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinError(8,mu21[16]->GetBinError(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinContent(9,mu21[25]->GetBinContent(1)*scaleFactor_mumu[8]);
  yieldPlot_mumu->SetBinError(9,mu21[25]->GetBinError(1)*scaleFactor_mumu[8]);
  
  TH1D * yieldPlot_ee = new TH1D("yieldPlot_ee","",10,-0.25,9.75);
  yieldPlot_ee->SetBinContent(1,e[0]->GetBinContent(1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinError(1,e[0]->GetBinError(1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinContent(2,e[1]->GetBinContent(1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinError(2,e[1]->GetBinError(1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinContent(3,e[2]->GetBinContent(1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinError(3,e[2]->GetBinError(1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinContent(4,e[3]->GetBinContent(1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinError(4,e[3]->GetBinError(1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinContent(5,e[4]->GetBinContent(1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinError(5,e[4]->GetBinError(1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinContent(6,e[5]->GetBinContent(1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinError(6,e[5]->GetBinError(1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinContent(7,e[15]->GetBinContent(1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinError(7,e[15]->GetBinError(1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinContent(8,e[16]->GetBinContent(1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinError(8,e[16]->GetBinError(1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinContent(9,e[25]->GetBinContent(1)*scaleFactor_ee[8]);
  yieldPlot_ee->SetBinError(9,e[25]->GetBinError(1)*scaleFactor_ee[8]);

  TH1D * yieldCombo = new TH1D("yieldCombo","",10,0.25,10.25);
  for(int i = 0; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){//subtract 1 for the pp bin
    float mu = yieldPlot_mumu->GetBinContent(i);
    float muErr = yieldPlot_mumu->GetBinError(i);
    float ept = yieldPlot_ee->GetBinContent(i);
    float eErr = yieldPlot_ee->GetBinError(i);

    //calculate a weighted mean
    float point = mu/(muErr*muErr)+ept/(eErr*eErr);
    float norm = 1.0/(muErr*muErr) + 1.0/(eErr*eErr);

    yieldCombo->SetBinContent(i, point/norm ); 
    yieldCombo->SetBinError(i, TMath::Sqrt(1.0/norm));
  }

  TH1D * ATLAS = new TH1D("ATLAS","",1,9,10);
  ATLAS->SetBinContent(1,0.3745*TMath::Power(10,-6));
  ATLAS->SetBinError(1,0.0085*TMath::Power(10,-6));
  ATLAS->SetMarkerStyle(8);
  ATLAS->SetMarkerColor(kGreen+1);
  ATLAS->SetLineColor(kGreen+1);

  for(int i = 1; i<11; i++){
    yieldPlot_mumu->GetXaxis()->SetBinLabel(i, labels[i-1]);
    yieldPlot_mumu->GetXaxis()->ChangeLabel(i,45);
  }
  
  yieldPlot_mumu24->SetMarkerStyle(8);
  yieldPlot_mumu24->SetMarkerColor(kBlue+1);
  yieldPlot_mumu24->SetMarkerSize(1.3);
  yieldPlot_mumu24->SetLineColor(kBlue+1);
  yieldPlot_mumu24->SetLineWidth(2);
  
  yieldPlot_mumu->SetMarkerStyle(8);
  yieldPlot_mumu->SetMarkerColor(kViolet+1);
  yieldPlot_mumu->SetMarkerSize(1.3);
  yieldPlot_mumu->SetLineColor(kViolet+1);
  yieldPlot_mumu->SetLineWidth(2);
  yieldPlot_mumu->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z}");
  yieldPlot_mumu->GetYaxis()->SetRangeUser(0.1*TMath::Power(10,-6),0.5*TMath::Power(10,-6));

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
  yieldPlot_mumu24->Draw("same");
  yieldPlot_ee->Draw("same");
  yieldCombo->Draw("same");
  ATLAS->Draw("same"); 
 
  yieldPlot_mumu->Print("All");
  yieldPlot_ee->Print("All");
  yieldCombo->Print("All");  

  TLegend * leg = new TLegend(0.21,0.22,0.81,0.43);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb, p_{T}^{l} > 20 GeV","");
  leg->AddEntry(yieldPlot_mumu24,"#mu^{+}#mu^{-} |#eta^{l}| < 2.4","p");
  leg->AddEntry(yieldPlot_mumu,"#mu^{+}#mu^{-} |#eta^{l}| < 2.1","p");
  leg->AddEntry(yieldPlot_ee,"e^{+}e^{-} |#eta^{l}| < 2.1 ","p");
  leg->AddEntry(yieldCombo,"combined |#eta^{l}| < 2.1","p");
  leg->AddEntry(ATLAS,"ATLAS pp |#eta^{l}| < 2.5","p");

  leg->Draw("same");

  c1->SaveAs("plots/prettyPlots/yields_Pretty.png");
  c1->SaveAs("plots/prettyPlots/yields_Pretty.pdf");
  c1->SaveAs("plots/prettyPlots/yields_Pretty.C");


  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: massPeakPlots <Z2EE file> <Z2mumu21 file> <Z2mumu24 file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu21 = argv[2];
  std::string Zmumu24 = argv[3];
   
  plotMassPeaks(Zee, Zmumu21, Zmumu24);
  return 0; 
}

#include "include/ptCorrector.h"
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

void prettyPlots(std::string Zee, std::string Zmumu21, std::string Zmumu24){
  Settings s = Settings();

  ptCorrector ptCorr = ptCorrector("resources/Z2mumu_Efficiencies.root","resources/Z2ee_EfficiencyMC_0.root");

  gStyle->SetErrorX(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //CentralityTool c = CentralityTool();
  //const int nBins = c.getNCentBins();

  //mumu channel first
  TFile * mu24 = TFile::Open(Zmumu24.c_str(),"read");
  TH1D * y_mu24 = (TH1D*)mu24->Get("yOS_minusAll_0_90");
  TH1D * pt_mu24 = (TH1D*)mu24->Get("pTOS_minusAll_0_90");

  pt_mu24->Divide( ptCorr.correction[0] );
 
  TFile * mu21 = TFile::Open(Zmumu21.c_str(),"read");
  TH1D * y_mu21 = (TH1D*)mu21->Get("yOS_minusAll_0_90");
  TH1D * pt_mu21 = (TH1D*)mu21->Get("pTOS_minusAll_0_90");
  
  pt_mu21->Divide( ptCorr.correction[2] );
  
  TFile * e = TFile::Open(Zee.c_str(),"read");
  TH1D * y_e = (TH1D*)e->Get("yOS_minusAll_0_90");
  TH1D * pt_e = (TH1D*)e->Get("pTOS_minusAll_0_90");
  
  pt_e->Divide( ptCorr.correction[1] );


  //make muon plot
  y_mu24->Scale(s.netLumi/(s.muLumi * s.Nmb));
  y_mu21->Scale(s.netLumi/(s.muLumi * s.Nmb));
  y_e->Scale(s.netLumi/(s.eLumi * s.Nmb));
  y_e->GetYaxis()->SetTitle("#frac{1}{N_{MB}} #frac{dN_{Z}}{dy}");

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  
  y_mu24->SetMarkerColor(kBlue);
  y_mu24->SetMarkerStyle(21);
  y_mu24->SetLineColor(kBlue);
  
  y_mu21->SetMarkerColor(kViolet+1);
  y_mu21->SetMarkerStyle(8);
  y_mu21->SetLineColor(kViolet+1);
  
  y_e->SetMarkerColor(kRed);
  y_e->SetMarkerStyle(25);
  y_e->SetLineColor(kRed);
  y_e->GetXaxis()->CenterTitle();
  y_e->GetYaxis()->CenterTitle();
  y_e->SetStats(0);  
  y_e->GetXaxis()->SetRangeUser(-2.4,2.4);
  y_e->GetYaxis()->SetRangeUser(0,y_e->GetMaximum()*1.4);

  
  y_e->Draw("p");
  y_mu21->Draw("p same");
  y_mu24->Draw("p same");

  TLegend *ly = new TLegend(0.625,0.725,0.875,0.875);
  ly->SetBorderSize(0);
  ly->SetFillStyle(0);
  ly->AddEntry((TObject*)0,"PbPb (0-90%)","");
  ly->AddEntry(y_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
  ly->AddEntry(y_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
  ly->AddEntry(y_e,"Z #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
  ly->Draw("same");
  
  c1->SaveAs("plots/prettyPlots/rapidity_Pretty.png"); 
  c1->SaveAs("plots/prettyPlots/rapidity_Pretty.pdf"); 
  c1->SaveAs("plots/prettyPlots/rapidity_Pretty.C"); 
  
  //make pt plot
  pt_mu24->Scale(s.netLumi/(s.muLumi * s.Nmb));
  pt_mu21->Scale(s.netLumi/(s.muLumi * s.Nmb));
  pt_e->Scale(s.netLumi/(s.eLumi * s.Nmb));
  pt_e->GetYaxis()->SetTitle("#frac{1}{N_{MB}} #frac{dN_{Z}}{dp_{T}} (GeV^{-1})");
  
  pt_mu24->SetMarkerColor(kBlue);
  pt_mu24->SetMarkerStyle(21);
  pt_mu24->SetLineColor(kBlue);
  
  pt_mu21->SetMarkerColor(kViolet+1);
  pt_mu21->SetMarkerStyle(8);
  pt_mu21->SetLineColor(kViolet+1);
  
  pt_e->SetMarkerColor(kRed);
  pt_e->SetMarkerStyle(25);
  pt_e->SetLineColor(kRed);
    
  TH1D * dummy = new TH1D("dummy",";p_{T} (GeV);#frac{1}{N_{MB}} #frac{dN_{Z}}{dp_{T}} (GeV^{-1})",2,0.1,200);
  dummy->SetBinContent(1,pt_mu24->GetMaximum());
  dummy->SetBinContent(2,pt_e->GetMinimum());
  dummy->SetLineColor(kWhite);
  dummy->SetMarkerColor(kWhite);
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  dummy->SetStats(0);
  dummy->Draw();
  
  pt_e->Draw("p same");
  pt_mu21->Draw("p same");
  pt_mu24->Draw("p same");

  c1->SetLogy();
  c1->SetLogx();

  TLegend *lpt = new TLegend(0.325,0.325,0.575,0.475);
  lpt->SetBorderSize(0);
  lpt->SetFillStyle(0);
  lpt->AddEntry((TObject*)0,"PbPb (0-90%)","");
  lpt->AddEntry(pt_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
  lpt->AddEntry(pt_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
  lpt->AddEntry(pt_e,"Z #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
  lpt->Draw("same");
  
  c1->SaveAs("plots/prettyPlots/pt_Pretty.png"); 
  c1->SaveAs("plots/prettyPlots/pt_Pretty.pdf"); 
  c1->SaveAs("plots/prettyPlots/pt_Pretty.C"); 


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
   
  prettyPlots(Zee, Zmumu21, Zmumu24);
  return 0; 
}

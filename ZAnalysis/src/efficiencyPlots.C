#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


void effPlots(std::string Zmumu, std::string Zee){
  TH1::SetDefaultSumw2();
//  Settings s = Settings();

//  CentralityTool c = CentralityTool();
//  const int nBins = c.getNCentBins();


  TEfficiency * eff_pt[2];
  TEfficiency * eff_y[2];
  TEfficiency * eff_phi[2];
  TEfficiency * eff_cent[2];
  TEfficiency * eff_2D[2];

  TGraphAsymmErrors * g_pt[2];
  TGraphAsymmErrors * g_y[2];
  TGraphAsymmErrors * g_phi[2];
  TGraphAsymmErrors * g_cent[2];

  TH2 * efficiency2D[2];


  TFile * mu = TFile::Open(Zmumu.c_str(), "read");
  TFile * ee = TFile::Open(Zee.c_str(), "read"); 

  eff_pt[0] = (TEfficiency*) mu->Get("eff_pt_0_90");
  eff_y[0] = (TEfficiency*) mu->Get("eff_y_0_90");
  eff_phi[0] = (TEfficiency*) mu->Get("eff_phi_0_90");
  eff_cent[0] = (TEfficiency*) mu->Get("eff_cent_0_100");
  eff_pt[1] = (TEfficiency*) ee->Get("eff_pt_0_90");
  eff_y[1] = (TEfficiency*) ee->Get("eff_y_0_90");
  eff_phi[1] = (TEfficiency*) ee->Get("eff_phi_0_90");
  eff_cent[1] = (TEfficiency*) ee->Get("eff_cent_0_100");

  eff_2D[0] = (TEfficiency*) mu->Get("eff_0_90");
  eff_2D[1] = (TEfficiency*) ee->Get("eff_0_90");

  for(int i = 0; i<2; i++){
    g_pt[i] = eff_pt[i]->CreateGraph();
    g_y[i] = eff_y[i]->CreateGraph();
    g_phi[i] = eff_phi[i]->CreateGraph();
    g_cent[i] = eff_cent[i]->CreateGraph();

    g_pt[i]->GetYaxis()->SetRangeUser(0,1);
    g_y[i]->GetYaxis()->SetRangeUser(0,1);
    g_phi[i]->GetYaxis()->SetRangeUser(0,1);
    g_cent[i]->GetYaxis()->SetRangeUser(0,1);
    
    g_pt[i]->GetYaxis()->SetTitle("efficiency");
    g_y[i]->GetYaxis()->SetTitle("efficiency");
    g_phi[i]->GetYaxis()->SetTitle("efficiency");
    g_cent[i]->GetYaxis()->SetTitle("efficiency");
    
    g_pt[i]->GetXaxis()->SetTitle("p_{T}");
    g_y[i]->GetXaxis()->SetTitle("y");
    g_phi[i]->GetXaxis()->SetTitle("#phi");
    g_cent[i]->GetXaxis()->SetTitle("centrality");

    if(i==1){
      g_pt[i]->SetMarkerColor(kRed+1);
      g_y[i]->SetMarkerColor(kRed+1);
      g_phi[i]->SetMarkerColor(kRed+1);
      g_cent[i]->SetMarkerColor(kRed+1);
      g_pt[i]->SetLineColor(kRed+1);
      g_y[i]->SetLineColor(kRed+1);
      g_phi[i]->SetLineColor(kRed+1);
      g_cent[i]->SetLineColor(kRed+1);
    }
  }

  efficiency2D[0] = eff_2D[0]->CreateHistogram();
  efficiency2D[1] = eff_2D[1]->CreateHistogram();
  for(int i = 0; i<2; i++){
    efficiency2D[i]->GetXaxis()->SetTitle("y");
    efficiency2D[i]->GetYaxis()->SetTitle("p_{T}");
  }

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  
  TLegend * leg = new TLegend(0.5,0.23,0.7,0.35);
  leg->SetBorderSize(0);
  leg->AddEntry(g_pt[0],"#mu^{+}#mu^{-}","pl");
  leg->AddEntry(g_pt[1],"e^{+}e^{-}","pl");

  g_pt[0]->Draw("AP");
  g_pt[1]->Draw("P same");
  leg->Draw("same");
  c1->SaveAs("plots/efficiencies/EffVsPt.png");
  c1->SaveAs("plots/efficiencies/EffVsPt.pdf");
  c1->SaveAs("plots/efficiencies/EffVsPt.C");
  c1->SetLogx();
  g_pt[0]->Draw("AP");
  g_pt[1]->Draw("P same");
  leg->Draw("same");
  c1->SaveAs("plots/efficiencies/EffVsPt_log.png");
  c1->SaveAs("plots/efficiencies/EffVsPt_log.pdf");
  c1->SaveAs("plots/efficiencies/EffVsPt_log.C");
  c1->SetLogx(0);
  
  g_y[0]->Draw("AP");
  g_y[1]->Draw("P same");
  leg->Draw("same");
  c1->SaveAs("plots/efficiencies/EffVsY.png");
  c1->SaveAs("plots/efficiencies/EffVsY.pdf");
  c1->SaveAs("plots/efficiencies/EffVsY.C");
  
  g_phi[0]->Draw("AP");
  g_phi[1]->Draw("P same");
  leg->Draw("same");
  c1->SaveAs("plots/efficiencies/EffVsPhi.png");
  c1->SaveAs("plots/efficiencies/EffVsPhi.pdf");
  c1->SaveAs("plots/efficiencies/EffVsPhi.C");
  
  g_cent[0]->Draw("AP");
  g_cent[1]->Draw("P same");
  leg->Draw("same");
  c1->SaveAs("plots/efficiencies/EffVsCent.png");
  c1->SaveAs("plots/efficiencies/EffVsCent.pdf");
  c1->SaveAs("plots/efficiencies/EffVsCent.C");

  efficiency2D[0]->GetZaxis()->SetRangeUser(0,1.0);
  efficiency2D[0]->Draw("colz");
  c1->SaveAs("plots/efficiencies/Eff2D_mu.png");
  c1->SaveAs("plots/efficiencies/Eff2D_mu.pdf");
  c1->SaveAs("plots/efficiencies/Eff2D_mu.C");
  c1->SetLogy();
  c1->SaveAs("plots/efficiencies/Eff2D_mu_log.png");
  c1->SaveAs("plots/efficiencies/Eff2D_mu_log.pdf");
  c1->SaveAs("plots/efficiencies/Eff2D_mu_log.C");
  c1->SetLogy(0);
  
  efficiency2D[1]->GetZaxis()->SetRangeUser(0,1.0);
  efficiency2D[1]->Draw("colz");
  c1->SaveAs("plots/efficiencies/Eff2D_e.png");
  c1->SaveAs("plots/efficiencies/Eff2D_e.pdf");
  c1->SaveAs("plots/efficiencies/Eff2D_e.C");
  c1->SetLogy();
  c1->SaveAs("plots/efficiencies/Eff2D_e_log.png");
  c1->SaveAs("plots/efficiencies/Eff2D_e_log.pdf");
  c1->SaveAs("plots/efficiencies/Eff2D_e_log.C");
  c1->SetLogy(0);

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: effPlots <Z2ee file> < Z2mumu file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu = argv[2];
   
  effPlots(Zmumu, Zee);
  return 0; 
}

#include "TH2D.h"
#include "../include/Settings.h"
#include "TCanvas.h"
#include <string>
#include "TFile.h"
#include "TLine.h"

void makeEMPlots(std::string inFile, bool isMu, bool isPeripheral, std::string outFile){
  TFile * f = TFile::Open(inFile.c_str(),"read");
  TH2D * h;
  if(!isPeripheral) h  = (TH2D*)f->Get("candAcoVsPt_0_90");
  if(isPeripheral)  h = (TH2D*)f->Get("candAcoVsPt_70_90");

  h->GetXaxis()->SetTitle("A");
  h->GetXaxis()->SetLabelSize(0.02);
  h->GetYaxis()->SetTitle("p_{T}");
  h->GetXaxis()->SetRangeUser(0,0.005);
  h->GetYaxis()->SetRangeUser(0,5.0);

  TCanvas * c1 = new TCanvas("c1","",800,800);

  h->SetStats(0);
  h->Draw("colz");

  TLine * l = new TLine(0,0,0,0);
  l->SetLineColor(kRed);
  l->SetLineWidth(2);

  Settings s = Settings();

  if(isMu){
    l->DrawLine(0,s.minPtCutForPhotons,s.acoCutForPhotons,s.minPtCutForPhotons);
    l->DrawLine(s.acoCutForPhotons,0,s.acoCutForPhotons,s.minPtCutForPhotons);
  } else {
    l->DrawLine(0,s.minPtCutForPhotons_ELE,s.acoCutForPhotons_ELE,s.minPtCutForPhotons_ELE);
    l->DrawLine(s.acoCutForPhotons_ELE,0,s.acoCutForPhotons_ELE,s.minPtCutForPhotons_ELE);
  }
  
  l->SetLineColor(kGreen+1);

  if(isMu){
    l->DrawLine(0,s.minPtCutForPhotonsU,s.acoCutForPhotonsU,s.minPtCutForPhotonsU);
    l->DrawLine(s.acoCutForPhotonsU,0,s.acoCutForPhotonsU,s.minPtCutForPhotonsU);
    l->DrawLine(0,s.minPtCutForPhotonsD,s.acoCutForPhotonsD,s.minPtCutForPhotonsD);
    l->DrawLine(s.acoCutForPhotonsD,0,s.acoCutForPhotonsD,s.minPtCutForPhotonsD);
  } else {
    l->DrawLine(0,s.minPtCutForPhotonsU_ELE,s.acoCutForPhotonsU_ELE,s.minPtCutForPhotonsU_ELE);
    l->DrawLine(s.acoCutForPhotonsU_ELE,0,s.acoCutForPhotonsU_ELE,s.minPtCutForPhotonsU_ELE);
    l->DrawLine(0,s.minPtCutForPhotonsD_ELE,s.acoCutForPhotonsD_ELE,s.minPtCutForPhotonsD_ELE);
    l->DrawLine(s.acoCutForPhotonsD_ELE,0,s.acoCutForPhotonsD_ELE,s.minPtCutForPhotonsD_ELE);
  }

  c1->SaveAs(Form("plots/photonBackground/%s.pdf",outFile.c_str()));
  c1->SaveAs(Form("plots/photonBackground/%s.png",outFile.c_str()));
  c1->SaveAs(Form("plots/photonBackground/%s.C",outFile.c_str()));


}

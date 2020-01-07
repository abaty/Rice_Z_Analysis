#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <iostream>
#include "TColor.h"
#include "TLegend.h"

void pdfProcess(){

  TFile * in = TFile::Open("../resources/Z2ee_EfficiencyMC_0.root","read");

  TH1D * h[12];
  for(int i = 0; i<12; i++){
    h[i] = (TH1D*) in->Get(Form("pdfIDVsPt_%d",i));
    h[i]->SetLineColor(kBlack);
    h[i]->SetLineWidth(1);
  }

  TH1D * d = new TH1D("dummy",";p_{T};fraction of processes",1,0.5,200);
  d->GetYaxis()->SetRangeUser(0,1.5);
  d->SetStats(0);

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  THStack * s = new THStack("stack","");
  //s->GetXaxis()->SetTitle("p_{T}");
  //s->GetYaxis()->SetTitle("fraction of processes");
  c1->SetLogx();
  //s->GetYaxis()->SetRangeUser(0,1.5);

  h[0]->SetFillColor(kBlue);
  h[1]->SetFillColor(kAzure-4);
  h[2]->SetFillColor(kCyan);
  h[3]->SetFillColor(kGreen+3);
  h[4]->SetFillColor(kGreen);
  h[5]->SetFillColor(kGreen-2);
  h[6]->SetFillColor(kRed);
  h[7]->SetFillColor(kRed-7);
  h[8]->SetFillColor(kOrange+1);
  h[9]->SetFillColor(kViolet);
  h[10]->SetFillColor(kYellow+1);
  //h[11]->SetFillColor(kGray+2);

  h[10]->Add(h[11]);
 
  for(int i = 0; i<11; i++){
    s->Add(h[i]);
  }
  d->Draw();
  s->Draw("same HIST");

  TLegend * l1 = new TLegend(0.3,0.675,0.55,0.875);
  l1->SetBorderSize(0);
  l1->AddEntry(h[0],"u#bar{u} #rightarrow Z","f");
  l1->AddEntry(h[1],"d#bar{d} #rightarrow Z","f");
  l1->AddEntry(h[2],"Q#bar{Q} #rightarrow Z","f");
  l1->AddEntry(h[3],"u#bar{u} #rightarrow Z + j","f");
  l1->AddEntry(h[4],"d#bar{d} #rightarrow Z + j","f");
  l1->AddEntry(h[5],"Q#bar{Q} #rightarrow Z + j","f");

  TLegend * l2 = new TLegend(0.6,0.675,0.875,0.875);
  l2->SetBorderSize(0);
  l2->AddEntry(h[6],"g+u #rightarrow Z + j","f");
  l2->AddEntry(h[7],"g+d #rightarrow Z + j","f");
  l2->AddEntry(h[8],"g+other #rightarrow Z + j","f");
  l2->AddEntry(h[9],"gq #rightarrow Z + 2j","f");
  l2->AddEntry(h[10],"Other","f");
  //l2->AddEntry(h[11],"Other","f");

  l1->Draw("same");
  l2->Draw("same");
 
  c1->RedrawAxis(); 
  c1->SaveAs("../plots/pdfStudy/pdfProcess.png");
  c1->SaveAs("../plots/pdfStudy/pdfProcess.pdf");
}

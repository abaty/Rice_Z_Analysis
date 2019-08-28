#include "TH2D.h"
#include "TCanvas.h"
#include <string>
#include "TFile.h"
#include "TLine.h"

void formatData(TH1D* h){
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(8);
  h->SetLineColor(kBlack);
  h->SetStats(0);
}

void formatMC(TH1D* h){
  h->SetFillColor(kOrange);
  h->SetLineColor(kBlack);
}

void makeLeptonPlots(std::string inFileData, std::string inFileMC , bool isMu){

  TFile * data = TFile::Open(inFileData.c_str(),"read");
  TFile * mc   = TFile::Open(inFileMC.c_str(), "read");

  TH1D * lepPt[2];
  TH1D * lepEta[2];
  TH1D * lepPhi[2];

  lepPt[0] = (TH1D*)  data->Get("lepPt");
  lepEta[0] = (TH1D*) data->Get("lepEta");
  lepPhi[0] = (TH1D*) data->Get("lepPhi");
  lepPt[1] = (TH1D*)  mc->Get("lepPt");
  lepEta[1] = (TH1D*) mc->Get("lepEta");
  lepPhi[1] = (TH1D*) mc->Get("lepPhi");

  for(int i = 0; i<2; i++){
    float entries = lepPt[i]->Integral();
    lepPt[i]->Scale(1.0/entries);

    entries = lepEta[i]->Integral();
    lepEta[i]->Scale(1.0/entries);
  
    entries = lepPhi[i]->Integral();
    lepPhi[i]->Scale(1.0/entries);
  }

  formatData(lepPt[0]);
  formatData(lepEta[0]);
  formatData(lepPhi[0]);
  formatMC(lepPt[1]);
  formatMC(lepEta[1]);
  formatMC(lepPhi[1]);

  lepPt[0]->GetYaxis()->SetTitle("Normalized to Unity");
  lepEta[0]->GetYaxis()->SetTitle("Normalized to Unity");
  lepPhi[0]->GetYaxis()->SetTitle("Normalized to Unity");

  TCanvas * c1 = new TCanvas("c1","",800,800);
  c1->SetLeftMargin(0.15);
  lepEta[0]->GetXaxis()->SetTitle("#eta_{l}");
  lepEta[0]->GetXaxis()->SetRangeUser(-2.4,2.4);
  lepEta[0]->Draw("p");
  lepEta[1]->Draw("h same"); 
  lepEta[0]->Draw("p same");

  TLegend * l = new TLegend(0.7,0.72,0.89,0.89);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(lepEta[0],"Z Data","p");
  l->AddEntry(lepEta[1],"Z MC","f");
  if(isMu) l->AddEntry((TObject*)0,"Z #rightarrow #mu^{+}#mu^{-}","");
  else l->AddEntry((TObject*)0,"Z #rightarrow e^{+}e^{-}","");
  
  l->Draw("same");

  c1->RedrawAxis();
  c1->SaveAs(Form("plots/leptonPlots/eta_isMu%d.png",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/eta_isMu%d.pdf",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/eta_isMu%d.C",(int)isMu));

  lepPhi[0]->GetXaxis()->SetTitle("#phi_{l}");
  lepPhi[0]->GetXaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
  lepPhi[0]->GetYaxis()->SetRangeUser(0,0.08);
  lepPhi[0]->Draw("p");
  lepPhi[1]->Draw("h same"); 
  lepPhi[0]->Draw("p same");

  l->Draw("same");
  c1->RedrawAxis();
  c1->SaveAs(Form("plots/leptonPlots/phi_isMu%d.png",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/phi_isMu%d.pdf",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/phi_isMu%d.C",(int)isMu));
  
  lepPt[0]->GetXaxis()->SetTitle("p_{T}^{l}");
  lepPt[0]->GetXaxis()->SetRangeUser(20,200);
  lepPt[0]->Draw("p");
  lepPt[1]->Draw("h same"); 
  lepPt[0]->Draw("p same");

  l->Draw("same");
  c1->RedrawAxis();
  c1->SetLogx();
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d.png",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d.pdf",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d.C",(int)isMu));
  c1->SetLogy();
  lepPt[0]->GetYaxis()->SetRangeUser(0.0001,1.0);
  lepPt[0]->Draw("p");
  lepPt[1]->Draw("h same"); 
  lepPt[0]->Draw("p same");
  l->Draw("same");
  c1->RedrawAxis();
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d_log.png",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d_log.pdf",(int)isMu));
  c1->SaveAs(Form("plots/leptonPlots/pt_isMu%d_log.C",(int)isMu));
}

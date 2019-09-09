#include "TCanvas.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <string>
#include "TFile.h"

void acceptancePlots(std::string inputFile, bool isMu, bool is21){

  TFile * input = TFile::Open(inputFile.c_str(),"read");

  TH1D * acceptance_y;
  TH1D * acceptance_pt;
  TH1D * variations_y[150];
  TH1D * variations_pt[150];  

  if(is21){
    acceptance_y =  (TH1D*) input->Get("accept21_y_ratio_0");
    acceptance_pt = (TH1D*) input->Get("accept21_pt_ratio_0");
  }else{
    acceptance_y =  (TH1D*) input->Get("accept24_y_ratio_0");
    acceptance_pt = (TH1D*) input->Get("accept24_pt_ratio_0");
  }

  for(int i = 0; i<105; i++){
    if(i==0 || (i>1 && i<9)) continue;
    if(is21){
      variations_y[i] = (TH1D*) input->Get(Form("accept21_nPDFVariation_y_%d",i));
      variations_pt[i] = (TH1D*) input->Get(Form("accept21_nPDFVariation_pt_%d",i));
    }
    else{
      variations_y[i] = (TH1D*) input->Get(Form("accept24_nPDFVariation_y_%d",i));
      variations_pt[i] = (TH1D*) input->Get(Form("accept24_nPDFVariation_pt_%d",i));
    }

    if(i==1) variations_y[i]->SetLineColor(kBlue);
    if(i==1) variations_pt[i]->SetLineColor(kBlue);
    if(i>1) variations_y[i]->SetLineColor(kGray);
    if(i>1) variations_pt[i]->SetLineColor(kGray);
    
  }

  TCanvas * cA = new TCanvas("cA","cA",800,800);
  acceptance_y->SetMarkerStyle(8);
  acceptance_y->SetLineColor(kBlack);
  acceptance_y->SetMarkerColor(kBlack);
  acceptance_y->GetYaxis()->SetRangeUser(0,1.0);
  acceptance_y->GetXaxis()->SetTitle("y_{Z}");
  acceptance_y->GetYaxis()->SetTitle("acceptance");
  acceptance_y->SetStats(0);
  acceptance_y->Draw("p");

  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_y_isMu%d.png", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_y_isMu%d.pdf", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_y_isMu%d.C", is21? 21 : 24, isMu? 1 : 0));

  cA->SetLeftMargin(0.2);
  for(int i = 0; i<105; i++){
    if(i==0 || (i>1 && i<9)) continue;
    if(i==1){
      variations_y[i]->GetXaxis()->SetTitle("y");
      variations_y[i]->SetStats(0);
      variations_y[i]->GetYaxis()->SetTitle("Ratio with Nominal Result");
      variations_y[i]->GetYaxis()->SetRangeUser(0.96,1.04);
      variations_y[i]->Draw("HIST l");
    }
    else variations_y[i]->Draw("HIST l same");
  }
  variations_y[1]->Draw("HIST l same");
  TLegend * l = new TLegend(0.6,0.6,0.9,0.9);
  l->AddEntry(variations_y[1],"nCTEQ15","l");
  l->AddEntry(variations_y[10],"EPPS16 copies","l");
  l->SetBorderSize(0);
  l->Draw("same");
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_y_isMu%d.png", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_y_isMu%d.pdf", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_y_isMu%d.C", is21? 21 : 24, isMu? 1 : 0));

  TH1D * dummy = new TH1D("dummy",";p_{T};acceptance",1,0.1,200);
  dummy->GetYaxis()->SetRangeUser(0,1.0);
  dummy->SetStats(0);
  dummy->Draw("");

  acceptance_pt->SetMarkerStyle(8);
  acceptance_pt->SetLineColor(kBlack);
  acceptance_pt->SetMarkerColor(kBlack);
  //acceptance_pt->GetXaxis()->SetTitle("p_{T}^{Z}");
  //acceptance_pt->GetYaxis()->SetTitle("acceptance");
  acceptance_pt->Draw("p same");
  cA->SetLogx();
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_pt_isMu%d.png", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_pt_isMu%d.pdf", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_pt_isMu%d.C", is21? 21 : 24, isMu? 1 : 0));
 
  dummy->GetYaxis()->SetTitle("Ratio with Nominal Result");
  dummy->GetYaxis()->SetRangeUser(0.96,1.04);
  dummy->SetStats(0);
  dummy->Draw(""); 
  cA->SetLeftMargin(0.2);
  for(int i = 0; i<105; i++){
    if(i==0 || (i>1 && i<9)) continue;
    if(i==1){
      variations_pt[i]->GetXaxis()->SetTitle("p_{T}");
      variations_pt[i]->SetStats(0);
      variations_pt[i]->GetYaxis()->SetTitle("Ratio with Nominal Result");
      variations_pt[i]->GetYaxis()->SetRangeUser(0.96,1.04);
      variations_pt[i]->Draw("HIST l same");
    }
    else variations_pt[i]->Draw("HIST l same");
  }
  variations_pt[1]->Draw("HIST l same");
  TLegend * l2 = new TLegend(0.6,0.6,0.9,0.9);
  l2->AddEntry(variations_pt[1],"nCTEQ15","l");
  l2->AddEntry(variations_pt[10],"EPPS16 copies","l");
  l2->SetBorderSize(0);
  l2->Draw("same");
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_pt_isMu%d.png", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_pt_isMu%d.pdf", is21? 21 : 24, isMu? 1 : 0));
  cA->SaveAs(Form("plots/acceptancePlots/acceptance%d_Syst_pt_isMu%d.C", is21? 21 : 24, isMu? 1 : 0));

}

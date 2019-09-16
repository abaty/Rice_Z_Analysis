#include "TFile.h"
#include "TH1D.h"
#include <string>

void makePlots(std::string inputFile, bool isMu){

  TFile * f = TFile::Open(inputFile.c_str(),"read");
  TH1D * eff = (TH1D*) f->Get("recoEff_pt_noPtWeight_0_90");
  TH1D * effWeighted = (TH1D*) f->Get("recoEff_pt_0_90");
  
  effWeighted->Divide(eff);

  effWeighted->GetYaxis()->SetRangeUser(0.95,1.05);
  effWeighted->GetYaxis()->SetTitle("Weighted/Unweighted efficiency");
  effWeighted->GetXaxis()->SetTitle("p_{T}");
  
  TCanvas * c1 = new TCanvas("c","c",800,800);
  effWeighted->SetStats(0);
  effWeighted->Draw();

  c1->SaveAs(Form("plots/ptSpectrumReweight/efficiencyReweightingComparison_isMu%d.png",(int) isMu));
  c1->SaveAs(Form("plots/ptSpectrumReweight/efficiencyReweightingComparison_isMu%d.pdf",(int) isMu));
  c1->SaveAs(Form("plots/ptSpectrumReweight/efficiencyReweightingComparison_isMu%d.C",(int) isMu));


  TH1D * anum = (TH1D*) f->Get(Form("accept%d_pt_noPtWeight_pass_0",isMu?24:21));
  TH1D * adenom = (TH1D*) f->Get(Form("accept%d_pt_noPtWeight_net_0",isMu?24:21));
  TH1D * aWeighted = (TH1D*) f->Get(Form("accept%d_pt_ratio_0",isMu?24:21));

  aWeighted->Divide(anum);
  aWeighted->Multiply(adenom);

  aWeighted->GetYaxis()->SetRangeUser(0.95,1.05);
  aWeighted->GetYaxis()->SetTitle("Wieghted/Unweighted acceptance");
  aWeighted->GetYaxis()->SetTitle("p_{T}");
  
  aWeighted->SetStats(0);
  aWeighted->Draw();


  c1->SaveAs(Form("plots/ptSpectrumReweight/acceptanceReweightingComparison_isMu%d.png",(int) isMu));
  c1->SaveAs(Form("plots/ptSpectrumReweight/acceptanceReweightingComparison_isMu%d.pdf",(int) isMu));
  c1->SaveAs(Form("plots/ptSpectrumReweight/acceptanceReweightingComparison_isMu%d.C",(int) isMu));



}

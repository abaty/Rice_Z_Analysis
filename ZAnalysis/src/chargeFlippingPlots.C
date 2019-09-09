#include <string>
#include "TH1D.h"
#include "TFile.h"


void chargeFlippingPlots(std::string inputData, std::string inputMC){

  TFile * data = TFile::Open(inputData.c_str(),"read");

  TH1D * sigData[3];
  sigData[0] = (TH1D*) data->Get("massPeakOS_withEff_0_90");
  sigData[1] = (TH1D*) data->Get("yOS_withEff_0_90");
  sigData[2] = (TH1D*) data->Get("pTOS_withEff_0_90");
  TH1D * ssData[3];
  ssData[0] = (TH1D*) data->Get("massPeakSS_withEff_0_90");
  ssData[1] = (TH1D*) data->Get("ySS_withEff_0_90");
  ssData[2] = (TH1D*) data->Get("pTSS_withEff_0_90");
  //TH1D * sigData = (TH1D*) data->Get("massPeakOS_minusAll_0_90");
  //TH1D * ssData = (TH1D*) data->Get("massPeakSS_ChargeFlipCorrected_0_90");


  TFile * mc = TFile::Open(inputMC.c_str(),"read");
  TH1D * sigMC[3];
  sigMC[0] = (TH1D*) mc->Get("massPeakOS_withEff_0_90");
  sigMC[1] = (TH1D*) mc->Get("yOS_withEff_0_90");
  sigMC[2] = (TH1D*) mc->Get("pTOS_withEff_0_90");
  TH1D * ssMC[3];
  ssMC[0] = (TH1D*) mc->Get("massPeakSS_withEff_0_90");
  ssMC[1] = (TH1D*) mc->Get("ySS_withEff_0_90");
  ssMC[2] = (TH1D*) mc->Get("pTSS_withEff_0_90");

  for(int i = 0; i<3; i++){
    ssData[i]->Divide(sigData[i]);
    ssMC[i]->Divide(sigMC[i]);
  }

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  for(int i = 0; i<3; i++){
    ssData[i]->SetMarkerStyle(8);
    ssData[i]->SetLineColor(kBlack);
    ssMC[i]->SetLineColor(kRed);

    ssData[i]->GetYaxis()->SetRangeUser(0,0.1);
    ssData[i]->GetYaxis()->SetTitle("Same sign / Opposite sign");

    TH1D * dummy;
    if(i==2){
      dummy = new TH1D("dummy",";m_{ee};Same sign / Opposite sign",1,0.1,200);
      dummy->GetYaxis()->SetRangeUser(0,0.1);
      dummy->SetStats(0);
      dummy->Draw();
    }
    

    ssData[i]->SetStats(0);
    if(i==2) c1->SetLogx();
    if(i!=2) ssData[i]->Draw("p");
    if(i==2) ssData[i]->Draw("p same");

    ssMC[i]->Draw("same h E0");

    TLegend * l = new TLegend(0.5,0.7,0.85,0.85);
    l->AddEntry(ssData[0],"Data (swapping + SS backgrounds)","p");
    l->AddEntry(ssMC[0],"MC (swapping only)","l");
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw("same");

    if(i==0)c1->SaveAs("plots/chargeSwapping/massPeak_chargeFlip.png");
    if(i==1)c1->SaveAs("plots/chargeSwapping/y_chargeFlip.png");
    if(i==2)c1->SaveAs("plots/chargeSwapping/pT_chargeFlip.png");
    
    if(i==0)c1->SaveAs("plots/chargeSwapping/massPeak_chargeFlip.pdf");
    if(i==1)c1->SaveAs("plots/chargeSwapping/y_chargeFlip.pdf");
    if(i==2)c1->SaveAs("plots/chargeSwapping/pT_chargeFlip.pdf");
    delete l;
  }
}

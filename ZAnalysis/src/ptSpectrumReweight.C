#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include <vector>
#include <string>

void normalize(TH1D* h){
  float entries = 0;
  entries = h->Integral();
  h->Scale(1.0/entries);
}

void ptSpectrumReweight(std::string eeData, std::string eeMC, std::string mumuData, std::string mumuMC){

  TFile * eD = TFile::Open(eeData.c_str(),"read");
  TFile * eMC = TFile::Open(eeMC.c_str(),"read");
  TFile * muD = TFile::Open(mumuData.c_str(),"read");
  TFile * muMC = TFile::Open(mumuMC.c_str(),"read");

  TH1D * eD_h = (TH1D*) eD->Get("candPt_withSFwithEvtWeight_0_90");
  TH1D * eMC_h = (TH1D*) eMC->Get("candPt_withSFwithEvtWeight_0_90");
  TH1D * muD_h = (TH1D*) muD->Get("candPt_withSFwithEvtWeight_0_90");
  TH1D * muMC_h = (TH1D*) muMC->Get("candPt_withSFwithEvtWeight_0_90");

  normalize(eD_h);
  normalize(eMC_h);
  normalize(muD_h);
  normalize(muMC_h);

  eD_h->Divide(eMC_h);
  muD_h->Divide(muMC_h);

  TCanvas * c1_e = new TCanvas("c1_e","",800,800);
  //c1_e->SetLogx();
  TH1D * dummy = new TH1D("dummy",";p_{T};Data/MC",1,0.1,200);
  dummy->GetYaxis()->SetRangeUser(0,2);
  dummy->Draw();
  dummy->SetStats(0);

  eD_h->SetMarkerStyle(8);
  eD_h->SetLineColor(kBlack);
  eD_h->SetLineColor(kBlack);
  eD_h->SetStats(0);
  eD_h->Draw("same A");
  
  TF1 * eleFit = new TF1("eleFit","[0]",70,200);
  eD_h->Fit("eleFit","EMRN0");
  double fitResult = eleFit->GetParameter(0);
  double fitErr = eleFit->GetParError(0);

  TH1D * ptWeightElectrons = new TH1D("ptWeightElectrons",";p_{T};",2000,0,200);
  TH1D * ptWeightElectronsU = new TH1D("ptWeightElectronsU",";p_{T};",2000,0,200);
  TH1D * ptWeightElectronsD = new TH1D("ptWeightElectronsD",";p_{T};",2000,0,200);
  for(int i = 1; i<2001; i++){
    float x = ptWeightElectrons->GetBinCenter(i);
    int xBin = eD_h->FindBin(x); 

    if(x<=80){
      float x1 = eD_h->GetBinCenter(xBin);
      float y1 = eD_h->GetBinContent(xBin);
      float y1U = eD_h->GetBinContent(xBin)+eD_h->GetBinError(xBin);
      float y1D = eD_h->GetBinContent(xBin)-eD_h->GetBinError(xBin);
      float x2 = 0, y2 = 0, y2U = 0, y2D = 0;
      float y3 = 0, y3U = 0, y3D = 0;
      if(x<=x1){
        x2 = eD_h->GetBinCenter(xBin-1);
        y2 = eD_h->GetBinContent(xBin-1);
        y2U = eD_h->GetBinContent(xBin-1) + eD_h->GetBinError(xBin-1);
        y2D = eD_h->GetBinContent(xBin-1) - eD_h->GetBinError(xBin-1);
        y3 = (y1-y2)/(x1-x2)*(x-x2) + y2;
        y3U = (y1U-y2U)/(x1-x2)*(x-x2) + y2U;
        y3D = (y1D-y2D)/(x1-x2)*(x-x2) + y2D; 
        if(x>=60){
          y1 = fitResult;
          x1 = 80;
          y3 = (y1-y2)/(x1-x2)*(x-x2) + y2;

          y1 = fitResult+fitErr; 
          y3U = (y1-y2U)/(x1-x2)*(x-x2) + y2U;
          
          y1 = fitResult-fitErr; 
          y3D = (y1-y2D)/(x1-x2)*(x-x2) + y2D;
        }
      }else{
        x2 = eD_h->GetBinCenter(xBin+1);
        y2 = eD_h->GetBinContent(xBin+1);
        y2U = eD_h->GetBinContent(xBin+1) + eD_h->GetBinError(xBin+1);
        y2D = eD_h->GetBinContent(xBin+1) - eD_h->GetBinError(xBin+1);
        y3 = (y2-y1)/(x2-x1)*(x-x1) + y1; 
        y3U = (y1U-y2U)/(x1-x2)*(x-x2) + y2U;
        y3D = (y1D-y2D)/(x1-x2)*(x-x2) + y2D; 
        if(x>=60){
          y2 = fitResult;
          x2 = 80;
          y3 = (y2-y1)/(x2-x1)*(x-x1) + y1; 
          
          y2 = fitResult+fitErr; 
          y3U = (y2-y1U)/(x2-x1)*(x-x1) + y1U;
          
          y2 = fitResult-fitErr; 
          y3D = (y2-y1D)/(x2-x1)*(x-x1) + y1D;
        }
      }
      ptWeightElectrons->SetBinContent(i,y3);
      ptWeightElectronsU->SetBinContent(i,y3U);
      ptWeightElectronsD->SetBinContent(i,y3D);

    }else{
      ptWeightElectrons->SetBinContent(i, fitResult);
      ptWeightElectronsU->SetBinContent(i, fitResult+fitErr);
      ptWeightElectronsD->SetBinContent(i, fitResult-fitErr);
    }
  }
  //fix the first 5 points
  for(int i = 1; i<6; i++){
    float mirror = ptWeightElectrons->GetBinContent(11-i);
    ptWeightElectrons->SetBinContent(i,eD_h->GetBinContent(1) + (eD_h->GetBinContent(1)-mirror));
    mirror = ptWeightElectronsU->GetBinContent(11-i);
    ptWeightElectronsU->SetBinContent(i,eD_h->GetBinContent(1)+eD_h->GetBinError(1) + (eD_h->GetBinContent(1)+eD_h->GetBinError(1)-mirror));
    mirror = ptWeightElectronsD->GetBinContent(11-i);
    ptWeightElectronsD->SetBinContent(i,eD_h->GetBinContent(1)-eD_h->GetBinError(1) + (eD_h->GetBinContent(1)-eD_h->GetBinError(1)-mirror));
  }

  ptWeightElectrons->SetLineColor(kRed);
  ptWeightElectronsU->SetLineColor(kBlue);
  ptWeightElectronsD->SetLineColor(kBlue);
  ptWeightElectrons->Draw("L same");
  ptWeightElectronsU->Draw("L same");
  ptWeightElectronsD->Draw("L same");

  c1_e->SaveAs("plots/ptSpectrumReweight/electron.png");
  c1_e->SaveAs("plots/ptSpectrumReweight/electron.pdf");
  c1_e->SaveAs("plots/ptSpectrumReweight/electron.C");
  
  TCanvas * c1_mu = new TCanvas("c1_mu","",800,800);
  //c1_mu->SetLogx();
  dummy->Draw();
  muD_h->SetMarkerStyle(8);
  muD_h->SetLineColor(kBlack);
  muD_h->SetLineColor(kBlack);
  muD_h->Draw("same");
  
  TF1 * muFit = new TF1("muFit","[0]",70,200);
  muD_h->Fit("muFit","EMRN0");
  fitResult = muFit->GetParameter(0);
  fitErr = muFit->GetParError(0);

  TH1D * ptWeightMuons = new TH1D("ptWeightMuons",";p_{T};",2000,0,200);
  TH1D * ptWeightMuonsU = new TH1D("ptWeightMuonsU",";p_{T};",2000,0,200);
  TH1D * ptWeightMuonsD = new TH1D("ptWeightMuonsD",";p_{T};",2000,0,200);
  for(int i = 1; i<2001; i++){
    float x = ptWeightMuons->GetBinCenter(i);
    int xBin = muD_h->FindBin(x); 

    if(x<=80){
      float x1 = muD_h->GetBinCenter(xBin);
      float y1 = muD_h->GetBinContent(xBin);
      float y1U = muD_h->GetBinContent(xBin)+muD_h->GetBinError(xBin);
      float y1D = muD_h->GetBinContent(xBin)-muD_h->GetBinError(xBin);
      float x2 = 0, y2 = 0, y2U = 0, y2D = 0;
      float y3 = 0, y3U = 0, y3D = 0;
      if(x<=x1){
        x2 = muD_h->GetBinCenter(xBin-1);
        y2 = muD_h->GetBinContent(xBin-1);
        y2U = muD_h->GetBinContent(xBin-1) + muD_h->GetBinError(xBin-1);
        y2D = muD_h->GetBinContent(xBin-1) - muD_h->GetBinError(xBin-1);
        y3 = (y1-y2)/(x1-x2)*(x-x2) + y2; 
        y3U = (y1U-y2U)/(x1-x2)*(x-x2) + y2U;
        y3D = (y1D-y2D)/(x1-x2)*(x-x2) + y2D; 
        if(x>=60){
          y1 = fitResult;
          x1 = 80;
          y3 = (y1-y2)/(x1-x2)*(x-x2) + y2; 
          
          y1 = fitResult+fitErr; 
          y3U = (y1-y2U)/(x1-x2)*(x-x2) + y2U;
          
          y1 = fitResult-fitErr; 
          y3D = (y1-y2D)/(x1-x2)*(x-x2) + y2D;
        }
      }else{
        x2 = muD_h->GetBinCenter(xBin+1);
        y2 = muD_h->GetBinContent(xBin+1);
        y2U = muD_h->GetBinContent(xBin+1) + muD_h->GetBinError(xBin+1);
        y2D = muD_h->GetBinContent(xBin+1) - muD_h->GetBinError(xBin+1);
        y3 = (y2-y1)/(x2-x1)*(x-x1) + y1; 
        y3U = (y2U-y1U)/(x2-x1)*(x-x1) + y1U;
        y3D = (y2D-y1D)/(x2-x1)*(x-x1) + y1D; 
        if(x>=60){
          y2 = fitResult;
          x2 = 80;
          y3 = (y1-y2)/(x1-x2)*(x-x2) + y2; 
          
          y2 = fitResult+fitErr; 
          y3U = (y2-y1U)/(x2-x1)*(x-x1) + y1U;
          
          y2 = fitResult-fitErr; 
          y3D = (y2-y1D)/(x2-x1)*(x-x1) + y1D;
        }
      }
      ptWeightMuons->SetBinContent(i,y3);
      ptWeightMuonsU->SetBinContent(i,y3U);
      ptWeightMuonsD->SetBinContent(i,y3D);

    }else{
      ptWeightMuons->SetBinContent(i, fitResult);
      ptWeightMuonsU->SetBinContent(i, fitResult+fitErr);
      ptWeightMuonsD->SetBinContent(i, fitResult-fitErr);
    }
  }
  //fix the first 5 points
  for(int i = 1; i<6; i++){
    float mirror = ptWeightMuons->GetBinContent(11-i);
    ptWeightMuons->SetBinContent(i,muD_h->GetBinContent(1) + (muD_h->GetBinContent(1)-mirror));
    mirror = ptWeightMuonsU->GetBinContent(11-i);
    ptWeightMuonsU->SetBinContent(i,muD_h->GetBinContent(1)+muD_h->GetBinError(1) + (muD_h->GetBinContent(1)+muD_h->GetBinError(1)-mirror));
    mirror = ptWeightMuonsD->GetBinContent(11-i);
    ptWeightMuonsD->SetBinContent(i,muD_h->GetBinContent(1)-muD_h->GetBinError(1) + (muD_h->GetBinContent(1)-muD_h->GetBinError(1)-mirror));
  }

  ptWeightMuons->Print("All");
  ptWeightMuons->SetLineColor(kRed);
  ptWeightMuonsU->SetLineColor(kBlue);
  ptWeightMuonsD->SetLineColor(kBlue);
  ptWeightMuons->Draw("L same");
  ptWeightMuonsU->Draw("L same");
  ptWeightMuonsD->Draw("L same");

  c1_mu->SaveAs("plots/ptSpectrumReweight/muon.png");
  c1_mu->SaveAs("plots/ptSpectrumReweight/muon.pdf");
  c1_mu->SaveAs("plots/ptSpectrumReweight/muon.C");

  TFile * output = TFile::Open("resources/ptSpectrumReweighting.root","recreate");
  ptWeightElectrons->Write();
  ptWeightElectronsU->Write();
  ptWeightElectronsD->Write();
  ptWeightMuons->Write();
  ptWeightMuonsU->Write();
  ptWeightMuonsD->Write();
}

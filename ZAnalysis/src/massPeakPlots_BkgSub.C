#include "include/RooUnfold/RooUnfoldResponse.h"
#include "include/RooUnfold/RooUnfoldInvert.h"
#include "include/RooUnfold/RooUnfold.h"
#include "include/CMS_lumi.C"
#include "include/HistNameHelper.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRandom3.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TLine.h"

void plotMassPeaks_BkgSub(std::string data_, std::string DY_, std::string ttbar_, std::string Wjet_, std::string responseFile_, bool isMu, std::string outTag){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  HistNameHelper h = HistNameHelper();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins][4][5]; 
  TH1D * massPeakSS[nBins][4][5]; 
  TH1D * massPeakSS_ChargeFlipCorrected[nBins][4][5]; 
  TH1D * massPeakOS_photons[nBins][4][5]; 
  TH1D * massPeakOS_minusSSAndPhoton[nBins][4][5];
  TH1D * massPeakOS_minusAll[nBins][4][5];

  TH1D * massPeakOS_DYsignal[nBins][4]; 
  TH1D * massPeakOS_DYphoton[nBins][4]; 
  TH1D * massPeakOS_DYtautau[nBins][4]; 
  TH1D * massPeakOS_DYsignalMinusPhoton[nBins][4];
  TH1D * massPeakOS_DYsignalMinusPhotonPlusSS[nBins][4];
  TH1D * massPeakOS_DYsignalMinusPhotonPlusBkg[nBins][4];
  TH1D * massPeakSS_DY[nBins][4];

  TH1D * massPeakOS_Wjet[nBins][4]; 
  TH1D * massPeakOS_ttbar[nBins][4]; 
 
  TH1D * MC_nEvents_DY;
  TH1D * MC_nEvents_Wjet;
  TH1D * MC_nEvents_ttbar;

  TH1D * fraction_tau[nBins][4];
  TH1D * fraction_Wjet[nBins][4];
  TH1D * fraction_ttbar[nBins][4];

  TH1D * bkg_tau[nBins][4][5];
  TH1D * bkg_Wjet[nBins][4][5];
  TH1D * bkg_ttbar[nBins][4][5];
 
  TFile * data = TFile::Open(data_.c_str(),"read");
  TFile * DY = TFile::Open(DY_.c_str(),"read");
  TFile * ttbar = TFile::Open(ttbar_.c_str(),"read");
  TFile * Wjet = TFile::Open(Wjet_.c_str(),"read");
  
  MC_nEvents_DY = (TH1D*)DY->Get("nEvents");  
  MC_nEvents_Wjet = (TH1D*)Wjet->Get("nEvents");  
  MC_nEvents_ttbar= (TH1D*)ttbar->Get("nEvents");  

  TFile * out = TFile::Open(Form("backgroundSubtraction_%s_isMu%d.root",outTag.c_str(),(int)isMu),"recreate");
  
  for(int i = 0; i<nBins; i++){
  for(int j = 0; j<4; j++){
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      massPeakOS[i][j][k] = (TH1D*)data->Get(Form("%sOS_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakSS[i][j][k] = (TH1D*)data->Get(Form("%sSS_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_photons[i][j][k] = (TH1D*)data->Get(Form("%sOS_ptLT0p5acoLT0p001_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    }

    massPeakOS_DYsignal[i][j] = (TH1D*)DY->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignal[i][j]->SetName(Form("%sOS_DYsignal_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_DY[i][j] = (TH1D*)DY->Get(Form("%sSS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_DY[i][j]->SetName(Form("%sSS_DYsignal_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphoton[i][j] = (TH1D*)DY->Get(Form("%sOS_ptLT0p5acoLT0p001_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphoton[i][j]->SetName(Form("%sOS_DY_ptLT0p5acoLT0p001_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYtautau[i][j] = (TH1D*)DY->Get(Form("%sOS_TauTau_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS_Wjet[i][j] = (TH1D*)Wjet->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbar[i][j] = (TH1D*)ttbar->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    massPeakOS_Wjet[i][j]->SetName(Form("%sOS_Wjet_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbar[i][j]->SetName(Form("%sOS_ttbar_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
 
    //do normalizations for MC
    float DY_norm = MC_nEvents_DY->Integral();
    massPeakOS_DYsignal[i][j]->Scale(s.DY_XS/DY_norm);
    massPeakOS_DYphoton[i][j]->Scale(s.DY_XS/DY_norm);
    massPeakOS_DYtautau[i][j]->Scale(s.DY_XS/DY_norm);
    massPeakSS_DY[i][j]->Scale(s.DY_XS/DY_norm);

    float Wjet_norm = MC_nEvents_Wjet->Integral();
    massPeakOS_Wjet[i][j]->Scale(s.Wjet_XS/Wjet_norm);

    float ttbar_norm = MC_nEvents_ttbar->Integral();
    massPeakOS_ttbar[i][j]->Scale(s.ttbar_XS/ttbar_norm);


    //Remove photon background from signal (these get cut in data as well)
    massPeakOS_DYsignalMinusPhoton[i][j] = (TH1D*)massPeakOS_DYsignal[i][j]->Clone(Form("%sOS_DYsignalMinusPhoton_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhoton[i][j]->Add(massPeakOS_DYphoton[i][j],-1);
    //add back the same sign pairs that we lose from charge-flipping (note that I think we should not do this now)
    massPeakOS_DYsignalMinusPhotonPlusSS[i][j] = (TH1D*)massPeakOS_DYsignalMinusPhoton[i][j]->Clone(Form("%sOS_DYsignalMinusPhotonPlusSS_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    //massPeakOS_DYsignalMinusPhotonPlusSS[i][j]->Add(massPeakSS_DY[i][j],1);

    //add the other backgrounds
    massPeakOS_DYsignalMinusPhotonPlusBkg[i][j] = (TH1D*)massPeakOS_DYsignalMinusPhotonPlusSS[i][j]->Clone(Form("%sOS_DYsignalMinusPhotonPlusBkg_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhotonPlusBkg[i][j]->Add(massPeakOS_DYtautau[i][j]);
    massPeakOS_DYsignalMinusPhotonPlusBkg[i][j]->Add(massPeakOS_ttbar[i][j]);
    massPeakOS_DYsignalMinusPhotonPlusBkg[i][j]->Add(massPeakOS_Wjet[i][j]);

    //Calculate the fractions for tautau, Wjet, TTbar
    //tautau
    fraction_tau[i][j] = (TH1D*)massPeakOS_DYtautau[i][j]->Clone(Form("%sFraction_tau_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_tau[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkg[i][j] );

    //ttbar
    fraction_ttbar[i][j] = (TH1D*)massPeakOS_ttbar[i][j]->Clone(Form("%sFraction_ttbar_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_ttbar[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkg[i][j] );
    
    //Wjet
    fraction_Wjet[i][j] = (TH1D*)massPeakOS_Wjet[i][j]->Clone(Form("%sFraction_Wjet_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_Wjet[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkg[i][j] );
   
    //remove entries not passing photon selection or SS from signal
    //data
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      massPeakOS_minusSSAndPhoton[i][j][k] = (TH1D*)massPeakOS[i][j][k]->Clone(Form("%sOS_minusSSAndPhoton%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_minusSSAndPhoton[i][j][k]->Add( massPeakOS_photons[i][j][k] , -1);
      
      //take charge flip into account
      massPeakSS_ChargeFlipCorrected[i][j][k] = (TH1D*) massPeakSS[i][j][k]->Clone(Form("%sSS_ChargeFlipCorrected%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
      //scan along bins looking for nonzero counts, and adjust based on MC estimation
      //commented lines are for adding the charge swaps back as signal (our MC efficiency assumes OS reco, so the same signs should not be added back in I guess...)
      for(int b = 0; b<massPeakSS_ChargeFlipCorrected[i][j][k]->GetSize(); b++){
        float SSYield = massPeakSS_ChargeFlipCorrected[i][j][k]->GetBinContent(b);
        if( SSYield == 0 ) continue;
        
        float estimate = ( massPeakSS_DY[i][j]->GetBinContent(b) / massPeakOS_DYsignalMinusPhotonPlusSS[i][j]->GetBinContent(b) ) * massPeakOS_minusSSAndPhoton[i][j][k]->GetBinContent(b);
        if( SSYield > estimate ){
          massPeakSS_ChargeFlipCorrected[i][j][k]->SetBinContent(b, SSYield-estimate);
         // massPeakOS_minusSSAndPhoton[i][j][k]->SetBinContent(b, massPeakOS_minusSSAndPhoton[i][j][k]->GetBinContent(b) + estimate );
        }
        else{
          massPeakSS_ChargeFlipCorrected[i][j][k]->SetBinContent(b, 0);
         // massPeakOS_minusSSAndPhoton[i][j][k]->SetBinContent(b, massPeakOS_minusSSAndPhoton[i][j][k]->GetBinContent(b) + SSYield );
        }
      }   
 
      //remove remaining SS background
      massPeakOS_minusSSAndPhoton[i][j][k]->Add( massPeakSS_ChargeFlipCorrected[i][j][k], -1);
    }    

    //normalize MC fractions to data with a ratio of counts between subtracted MC data
    //tau
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      bkg_tau[i][j][k] = (TH1D*) fraction_tau[i][j]->Clone(Form("%sBkg_tau%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_tau[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);
    
      //ttbar
      bkg_ttbar[i][j][k] = (TH1D*) fraction_ttbar[i][j]->Clone(Form("%sBkg_ttbar%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_ttbar[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);

      //Wjet
      bkg_Wjet[i][j][k] = (TH1D*) fraction_Wjet[i][j]->Clone(Form("%sBkg_Wjet%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_Wjet[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);

    //remove these backgrounds
      massPeakOS_minusAll[i][j][k] = (TH1D*) massPeakOS_minusSSAndPhoton[i][j][k]->Clone(Form("%sOS_minusAll%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_minusAll[i][j][k]->Add( bkg_tau[i][j][k], -1);
      massPeakOS_minusAll[i][j][k]->Add( bkg_ttbar[i][j][k], -1);
      massPeakOS_minusAll[i][j][k]->Add( bkg_Wjet[i][j][k], -1);
    }

    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      massPeakOS[i][j][k]->Write();
      massPeakSS[i][j][k]->Write();
      massPeakSS_ChargeFlipCorrected[i][j][k]->Write();
      massPeakOS_photons[i][j][k]->Write();
      massPeakOS_minusSSAndPhoton[i][j][k]->Write();
      massPeakOS_minusAll[i][j][k]->Write();
      bkg_Wjet[i][j][k]->Write(); 
      bkg_ttbar[i][j][k]->Write(); 
      bkg_tau[i][j][k]->Write(); 
    }
    massPeakOS_DYsignalMinusPhoton[i][j]->Write();
    fraction_Wjet[i][j]->Write(); 
    fraction_ttbar[i][j]->Write(); 
    fraction_tau[i][j]->Write(); 
    massPeakOS_DYsignal[i][j]->Write();
    massPeakOS_DYphoton[i][j]->Write();
    massPeakOS_DYtautau[i][j]->Write();
    massPeakOS_Wjet[i][j]->Write();
    massPeakOS_ttbar[i][j]->Write();

    //Draw the histogram Stack
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    //gStyle->SetErrorX(0);
    TCanvas * c1 = new TCanvas("c1","c1",800,900);
    c1->SetBorderSize(0);
    TPad * p1 = new TPad("p1","p1",0,0.2,1,1,0);
    TPad * p2 = new TPad("p2","p2",0,0,1,0.2,0);
    c1->SetLineWidth(0);
    p1->SetBottomMargin(0);
    p1->SetLeftMargin(0.15);
    p1->SetRightMargin(0.05);
    p1->SetTopMargin(0.07);
    p1->SetBorderSize(0);
    p1->Draw();
    p2->SetTopMargin(0);
    p2->SetLeftMargin(0.15);
    p2->SetRightMargin(0.05);
    p2->SetBottomMargin(0.3);
    p2->SetBorderSize(0);
    p2->Draw();
    p1->cd();
    
    if(j==0 && isMu) massPeakOS[i][j][0]->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
    if(j==0 && !isMu) massPeakOS[i][j][0]->GetXaxis()->SetTitle("m_{ee} (GeV)");
    if(j==1) massPeakOS[i][j][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
    if(j==2) massPeakOS[i][j][0]->GetXaxis()->SetTitle("y");
    if(j==3) massPeakOS[i][j][0]->GetXaxis()->SetTitle("yield");
    
    if(j==0 || j==3) massPeakOS[i][j][0]->GetYaxis()->SetTitle("Entries");
    if(j==1) massPeakOS[i][j][0]->GetYaxis()->SetTitle("#frac{dN_{Z}}{dp_{T}} (GeV^{-1})");
    if(j==2) massPeakOS[i][j][0]->GetYaxis()->SetTitle("#frac{dN_{Z}}{dy}");
    massPeakOS[i][j][0]->GetXaxis()->CenterTitle();
    massPeakOS[i][j][0]->GetYaxis()->CenterTitle();
    massPeakOS[i][j][0]->SetFillColor(kOrange+1);
    massPeakOS[i][j][0]->SetMarkerStyle(8);
    massPeakOS[i][j][0]->SetMarkerColor(kBlack);
    massPeakOS[i][j][0]->SetLineColor(kBlack);
    massPeakOS[i][j][0]->SetStats(0);

    TH1D * dummy;
    TH1D * dummyR;
    if(j==1){
      dummy = new TH1D("dummy",";p_{T} (GeV); #frac{dN_{Z}}{dp_{T}} (GeV^{-1})",2,0.1,200);
      dummyR = new TH1D("dummy",";p_{T} (GeV); #frac{Data}{MC}",2,0.1,200);
      dummy->SetBinContent(1,massPeakOS[i][j][0]->GetMaximum());
      dummy->SetBinContent(2,massPeakOS[i][j][0]->GetMinimum());
      dummy->SetLineColor(kWhite);
      dummy->GetYaxis()->CenterTitle();
      dummy->SetMarkerColor(kWhite);
      dummy->SetStats(0);
      dummyR->SetStats(0);
      dummy->Draw();
    }

    massPeakOS[i][j][0]->Draw("p same");

    bkg_ttbar[i][j][0]->SetFillColor(kGray);
    bkg_ttbar[i][j][0]->SetLineColor(kBlack);

    bkg_Wjet[i][j][0]->Add(bkg_ttbar[i][j][0]);
    bkg_Wjet[i][j][0]->SetFillColor(kViolet+1);
    bkg_Wjet[i][j][0]->SetLineColor(kBlack);
    
    bkg_tau[i][j][0]->Add(bkg_Wjet[i][j][0]);
    bkg_tau[i][j][0]->SetFillColor(kGreen+1);
    bkg_tau[i][j][0]->SetLineColor(kBlack);
    
    massPeakOS_photons[i][j][0]->Add(bkg_tau[i][j][0]);
    massPeakOS_photons[i][j][0]->SetFillColor(kRed+1);
    massPeakOS_photons[i][j][0]->SetLineColor(kBlack);
    
    massPeakSS_ChargeFlipCorrected[i][j][0]->Add(massPeakOS_photons[i][j][0]);
    massPeakSS_ChargeFlipCorrected[i][j][0]->SetFillColor(kBlue);
    massPeakSS_ChargeFlipCorrected[i][j][0]->SetLineColor(kBlack); 

    //normalize our signal to the number of data events minus the ones already accounted for
    float integralDYSignal = massPeakOS_DYsignalMinusPhoton[i][j]->Integral();
    float dataIntegral = massPeakOS[i][j][0]->Integral() - massPeakSS_ChargeFlipCorrected[i][j][0]->Integral();
    massPeakOS_DYsignalMinusPhoton[i][j]->Scale(dataIntegral/integralDYSignal);
    massPeakOS_DYsignalMinusPhoton[i][j]->Add(massPeakSS_ChargeFlipCorrected[i][j][0]);
    massPeakOS_DYsignalMinusPhoton[i][j]->SetFillColor(kOrange+1);
    massPeakOS_DYsignalMinusPhoton[i][j]->SetLineColor(kBlack);

    massPeakOS_DYsignalMinusPhoton[i][j]->Draw("HIST same");
    massPeakSS_ChargeFlipCorrected[i][j][0]->Draw("HIST same");
    massPeakOS_photons[i][j][0]->Draw("HIST same");
    bkg_tau[i][j][0]->Draw("HIST same");
    bkg_Wjet[i][j][0]->Draw("HIST same");
    bkg_ttbar[i][j][0]->Draw("HIST same");
    massPeakOS[i][j][0]->Draw("p same");

    TLegend *leg = new TLegend(0.65,0.65,0.92,0.885);
    leg->AddEntry(massPeakOS[i][j][0],Form("Data (%d-%d%%)",c.getCentBinLow(i),c.getCentBinHigh(i)),"p");
    if(isMu) leg->AddEntry(massPeakOS_DYsignalMinusPhoton[i][j],"Z #rightarrow #mu^{+}#mu^{-}","f");
    if(!isMu) leg->AddEntry(massPeakOS_DYsignalMinusPhoton[i][j],"Z #rightarrow e^{+}e^{-}","f");
    leg->AddEntry(massPeakSS_ChargeFlipCorrected[i][j][0],"Same Sign (QCD)","f");
    leg->AddEntry(massPeakOS_photons[i][j][0],"EM background","f");
    leg->AddEntry(bkg_tau[i][j][0],"Z #rightarrow #tau^{+}#tau^{-}","f");
    leg->AddEntry(bkg_Wjet[i][j][0],"W^{#pm} + X","f");
    leg->AddEntry(bkg_ttbar[i][j][0],"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    //ratio plot
    TH1D * tempRatio = (TH1D*)massPeakOS[i][j][0]->Clone("tempRatio");
    tempRatio->Divide(massPeakOS_DYsignalMinusPhoton[i][j]);
    tempRatio->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
    p2->cd();
    if(j==1){
      dummyR->GetYaxis()->SetTitleOffset(0.4);
      dummyR->GetYaxis()->SetRangeUser(0,1.99);
      dummyR->GetYaxis()->SetTitleSize(0.14);
      dummyR->GetXaxis()->SetTitleSize(0.14);
      dummyR->GetYaxis()->SetLabelSize(0.14);
      dummyR->GetXaxis()->SetLabelSize(0.14);
      dummyR->GetYaxis()->CenterTitle();
      dummyR->GetXaxis()->CenterTitle();
      dummyR->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
      dummyR->Draw();
    }
    tempRatio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    tempRatio->GetYaxis()->SetTitleOffset(0.4);
    tempRatio->GetYaxis()->SetRangeUser(0,1.99);
    tempRatio->GetYaxis()->SetTitleSize(0.14);
    tempRatio->GetXaxis()->SetTitleSize(0.14);
    tempRatio->GetYaxis()->SetLabelSize(0.14);
    tempRatio->GetXaxis()->SetLabelSize(0.14);
    tempRatio->Draw("same");
 
    TLine * line1;
    line1 = new TLine(tempRatio->GetXaxis()->GetBinLowEdge(1),1,tempRatio->GetXaxis()->GetBinUpEdge(tempRatio->GetSize()-2),1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");

    p1->cd();

    if(j==1){
      p2->SetLogx();
      tempRatio->GetXaxis()->SetRangeUser(0.1,199);
      p1->SetLogx();
      p1->SetLogy();
      dummy->GetXaxis()->SetRangeUser(0.1,199);
      dummy->GetYaxis()->SetRangeUser(0.01,200000);
    }

    p1->RedrawAxis();

    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));

    if(j==0) massPeakOS[i][j][0]->GetYaxis()->SetRangeUser(0.015, massPeakOS[i][j][0]->GetMaximum()*50);
    if(j==2) massPeakOS[i][j][0]->GetYaxis()->SetRangeUser(0.3, massPeakOS[i][j][0]->GetMaximum()*300);
    p1->SetLogy();
    CMS_lumi(p1,0,10,1.8);

    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
   
    if(j==1) delete dummy;
    if(j==1) delete dummyR;
    delete line1; 
    delete tempRatio;
    delete p1;
    delete p2;
    delete c1;
    delete leg;
  }
  }

 
  //forwarding some histograms into the output 
  TH1D * mcStatRelErr[nBins];
  TH1D * mcStatRelErr_pt;
  TH1D * mcStatRelErr_y;
  for(int i = 0; i<nBins; i++){
    mcStatRelErr[i] = (TH1D*) data->Get(Form("yieldOS_withEff_RelStatErr_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    out->cd();
    mcStatRelErr[i]->Write();
  } 
  mcStatRelErr_pt = (TH1D*)data->Get("pTOS_withEff_RelStatErr_0_90");
  mcStatRelErr_pt->Write();
  mcStatRelErr_y = (TH1D*)data->Get("yOS_withEff_RelStatErr_0_90");
  mcStatRelErr_y->Write();


  //UNFOLDING
  const int nToys = 1000;
  TFile * responseFile = TFile::Open(responseFile_.c_str(),"read");
  TH2 * unf_response = (TH2*) responseFile->Get("unfolding_response");
  TH2 * unf_responseU = (TH2*) responseFile->Get("unfolding_responseU");
  TH2 * unf_responseD = (TH2*) responseFile->Get("unfolding_responseD");
  TH2D * unf_responseHist = (TH2D*) responseFile->Get("unfolding_response");
  TH1D * unf_reco = (TH1D*) responseFile->Get("unfolding_recoPt");
  TH1D * unf_gen = (TH1D*) responseFile->Get("unfolding_genPt");

  TH2 * unf_responseLargeErrors[nToys];
  TRandom3 randGen = TRandom3();
  for(int j = 0; j<nToys; j++){
    unf_responseLargeErrors[j] = (TH2*) unf_response->Clone("unfolding_response");
    for(int i = 0; i<unf_responseHist->GetSize(); i++){
      float mean = unf_responseLargeErrors[j]->GetBinContent(i);
      float err = unf_responseLargeErrors[j]->GetBinError(i);
      float newPoint = randGen.Gaus(mean, err);
      unf_responseLargeErrors[j]->SetBinContent(i, newPoint);
    }
  }

  out->cd();
  RooUnfoldResponse nominalResponse = RooUnfoldResponse( 0, 0, unf_response);
  RooUnfoldResponse nominalResponseU = RooUnfoldResponse( 0, 0, unf_responseU);
  RooUnfoldResponse nominalResponseD = RooUnfoldResponse( 0, 0, unf_responseD);
  RooUnfoldResponse nominalResponseLargeErr[nToys];
  for(int i = 0; i<nToys; i++){
    nominalResponseLargeErr[i] = RooUnfoldResponse(0,0,unf_responseLargeErrors[i]);
  }

  h.undoDifferential( massPeakOS_minusAll[25][1][0] );
  RooUnfoldInvert unfoldObject = RooUnfoldInvert(&nominalResponse, massPeakOS_minusAll[25][1][0]);
  TH1D * unfoldedDist = (TH1D*) unfoldObject.Hreco();
  unfoldedDist->Print("All");
  RooUnfoldInvert unfoldObjectU = RooUnfoldInvert(&nominalResponseU, massPeakOS_minusAll[25][1][0]);
  TH1D * unfoldedDistU = (TH1D*) unfoldObjectU.Hreco();
  unfoldedDistU->Divide(unfoldedDist);
  unfoldedDistU->SetName("unfoldedDistU");
  unfoldedDistU->Print("All");
  RooUnfoldInvert unfoldObjectD = RooUnfoldInvert(&nominalResponseD, massPeakOS_minusAll[25][1][0]);
  TH1D * unfoldedDistD = (TH1D*) unfoldObjectD.Hreco();
  unfoldedDistD->Divide(unfoldedDist);
  unfoldedDistD->SetName("unfoldedDistD");
  unfoldedDistD->Print("All");
  TH1D * unfoldingUncert = (TH1D*)unfoldedDistU->Clone("unfoldingUncert_SpectrumShape");
  unfoldingUncert->Reset();
  for(int i = 0; i<unfoldingUncert->GetSize();i++){
    unfoldingUncert->SetBinContent(i, TMath::Max( TMath::Abs(unfoldedDistU->GetBinContent(i)-1), TMath::Abs(unfoldedDistD->GetBinContent(i)-1)));
  }

  //errors check
  TH1D * unfoldingUncert_Stats = (TH1D*)unfoldedDistU->Clone("unfoldingUncert_Stats");
  unfoldingUncert_Stats->Reset();
  TH1D * unfoldedDist_LargeErrors[nToys];
  RooUnfoldInvert unfoldObjectLargeErrors[nToys];
  for(int i = 0; i<nToys; i++){
    unfoldObjectLargeErrors[i] = RooUnfoldInvert(&nominalResponseLargeErr[i], massPeakOS_minusAll[25][1][0]);
    unfoldedDist_LargeErrors[i] = (TH1D*) unfoldObjectLargeErrors[i].Hreco();
    unfoldedDist_LargeErrors[i]->Divide(unfoldedDist);
  }
  unfoldedDist_LargeErrors[0]->Print("All");
  for(int i = 0; i<unfoldingUncert_Stats->GetSize();i++){
    std::vector< float > values;
    //float max = 0;
    for(int j = 0; j<nToys; j++) values.push_back( unfoldedDist_LargeErrors[j]->GetBinContent(i) );
    std::sort( values.begin(), values.end() );
    unfoldingUncert_Stats->SetBinContent(i, TMath::Max( TMath::Abs( values.at(((int)(0.683*nToys)))-1), TMath::Abs( values.at((int)((1-0.683)*nToys))-1) ) );
  }
  unfoldingUncert->Print("All");
  unfoldingUncert_Stats->Print("All");
  TH1D * unfoldingUncert_total = (TH1D*)unfoldingUncert->Clone("unfoldingUncert");
  h.addInQuadrature2( unfoldingUncert_total, unfoldingUncert_Stats);  
  unfoldingUncert_total->Print("All");  
   
  //smearing check
  TH1D * unfoldedDistWithSmearing = (TH1D*) unfoldedDist->Clone("unfoldedDistWithSmearing");
  unfoldedDistWithSmearing->Reset();
  for(int i = 1; i<unfoldedDist->GetSize()-1; i++){
    double rowSum = 0;
    for(int j = 1; j<unfoldedDist->GetSize()-1; j++) rowSum += unf_responseHist->GetBinContent(j, i);
    for(int j = 1; j<unfoldedDist->GetSize()-1; j++){
      double old = unfoldedDistWithSmearing->GetBinContent(j);
      unfoldedDistWithSmearing->SetBinContent(j, old + unfoldedDist->GetBinContent(i)*unf_responseHist->GetBinContent(j,i)/rowSum);
    } 
  }
  unfoldedDistWithSmearing->Print("All");
  unfoldedDistWithSmearing->Divide(massPeakOS_minusAll[25][1][0]);
  unfoldedDistWithSmearing->Print("All");

  h.makeDifferential( massPeakOS_minusAll[25][1][0] );
  massPeakOS_minusAll[25][1][0]->Print("All");

  unfoldedDist->SetName("unfoldedDist");
  h.makeDifferential(unfoldedDist);
  unfoldedDist->Print("All");

  TH1D * unfoldedRatio = (TH1D*) unfoldedDist->Clone("unfoldedRatio");
  unfoldedRatio->Divide(massPeakOS_minusAll[25][1][0]);
  unfoldedRatio->Print("All");
 
  //closure
  RooUnfoldInvert unfoldObjectClosure = RooUnfoldInvert(&nominalResponse, unf_reco);
  TH1D * unfoldedClosure = (TH1D*) unfoldObjectClosure.Hreco();
  unfoldedClosure->SetName("unfoldedClosure");
  unfoldedClosure->Divide(unf_gen); 
  unfoldedClosure->Print("All");

  
  TCanvas * respC = new TCanvas("respC","respC",800,800);
  respC->SetLeftMargin(0.2);
  unf_responseHist->GetXaxis()->SetTitle("p_{T}^{reco}");
  unf_responseHist->GetYaxis()->SetTitle("p_{T}^{gen}");
  unf_responseHist->SetStats(0);
  unf_responseHist->Draw("colz");
  respC->SetLogz();
  respC->SaveAs(Form("plots/Unfolding/response_isMu%d_%s.png",(int)isMu,outTag.c_str()));
  respC->SaveAs(Form("plots/Unfolding/response_isMu%d_%s.pdf",(int)isMu,outTag.c_str()));

  TH1D * dummy4 = new TH1D("dummy4",";p_{T};Unfolded/Raw",1,0.1,200);
  dummy4->Draw();
  dummy4->SetStats(0);
  respC->SetLogx();
  dummy4->GetYaxis()->SetRangeUser(0.8,1.2);
  unfoldedRatio->SetLineColor(kBlack);
  unfoldedRatio->SetMarkerStyle(8);
  unfoldedRatio->SetMarkerColor(kBlack);
  unfoldedRatio->Draw("same");
  respC->SaveAs(Form("plots/Unfolding/UnfoldingRatio_isMu%d_%s.png",(int)isMu,outTag.c_str()));
  respC->SaveAs(Form("plots/Unfolding/UnfoldingRatio_isMu%d_%s.pdf",(int)isMu,outTag.c_str()));
  
  dummy4->GetYaxis()->SetTitle("varied/nominal");
  dummy4->Draw();
  for(int i = 0; i<nToys; i++){
    unfoldedDist_LargeErrors[i]->SetLineColor(kGray);
    unfoldedDist_LargeErrors[i]->Draw("same HIST L");
  }
  respC->SaveAs(Form("plots/Unfolding/UnfoldingStatVariations_isMu%d_%s.png",(int)isMu,outTag.c_str()));
  respC->SaveAs(Form("plots/Unfolding/UnfoldingStatVariations_isMu%d_%s.pdf",(int)isMu,outTag.c_str()));

  dummy4->GetYaxis()->SetTitle("(unfolded+smeared) / raw");
  dummy4->Draw();
  unfoldedDistWithSmearing->Draw("same");
  respC->SaveAs(Form("plots/Unfolding/UnfoldingSmeared_isMu%d_%s.png",(int)isMu,outTag.c_str()));
  respC->SaveAs(Form("plots/Unfolding/UnfoldingSmeared_isMu%d_%s.pdf",(int)isMu,outTag.c_str()));
 
  dummy4->GetYaxis()->SetTitle("unfolded / gen");
  dummy4->Draw();
  unfoldedClosure->Draw("same");
  respC->SaveAs(Form("plots/Unfolding/UnfoldingClosure_isMu%d_%s.png",(int)isMu,outTag.c_str()));
  respC->SaveAs(Form("plots/Unfolding/UnfoldingClosure_isMu%d_%s.pdf",(int)isMu,outTag.c_str()));

 
  unfoldedDist->Write();
  unfoldedRatio->Write();
  unfoldedClosure->Write();
  unfoldedDistWithSmearing->Write();
  unfoldingUncert->Write();
  unfoldingUncert_Stats->Write();
  unfoldingUncert_total->Write();

  out->Close();
  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 8)
  {
    std::cout << "Usage: massPeakPlots_BkgSub <Data file> <DY File> <TTbar File> <WJet File> <ResponseMatrix file> <isMu> <outTag>" << std::endl;
    return 1;
  }  

  std::string data = argv[1];
  std::string DY = argv[2];
  std::string ttbar = argv[3];
  std::string Wjet = argv[4];
  std::string responseMatrix = argv[5];   
  bool isMu = (bool)std::atoi(argv[6]);
  std::string outTag = argv[7];

  plotMassPeaks_BkgSub(data, DY, ttbar, Wjet, responseMatrix, isMu, outTag);
  return 0; 
}

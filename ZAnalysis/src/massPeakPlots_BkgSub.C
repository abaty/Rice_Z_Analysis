#include "include/HistNameHelper.h"
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

void plotMassPeaks_BkgSub(std::string data_, std::string DY_, std::string ttbar_, std::string Wjet_, bool isMu, std::string outTag){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  HistNameHelper h = HistNameHelper();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins][4]; 
  TH1D * massPeakSS[nBins][4]; 
  TH1D * massPeakOS_photons[nBins][4]; 
  TH1D * massPeakOS_minusSSAndPhoton[nBins][4];
  TH1D * massPeakOS_minusAll[nBins][4];

  TH1D * massPeakOS_DYsignal[nBins][4]; 
  TH1D * massPeakOS_DYphoton[nBins][4]; 
  TH1D * massPeakOS_DYtautau[nBins][4]; 
  TH1D * massPeakOS_DYsignalMinusPhoton[nBins][4];
  TH1D * massPeakOS_DYsignalMinusPhotonPlusBkg[nBins][4];

  TH1D * massPeakOS_Wjet[nBins][4]; 
  TH1D * massPeakOS_ttbar[nBins][4]; 
 
  TH1D * MC_nEvents_DY;
  TH1D * MC_nEvents_Wjet;
  TH1D * MC_nEvents_ttbar;

  TH1D * fraction_tau[nBins][4];
  TH1D * fraction_Wjet[nBins][4];
  TH1D * fraction_ttbar[nBins][4];

  TH1D * bkg_tau[nBins][4];
  TH1D * bkg_Wjet[nBins][4];
  TH1D * bkg_ttbar[nBins][4];
 
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
    massPeakOS[i][j] = (TH1D*)data->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS[i][j] = (TH1D*)data->Get(Form("%sSS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_photons[i][j] = (TH1D*)data->Get(Form("%sOS_ptLT0p5acoLT0p001_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS_DYsignal[i][j] = (TH1D*)DY->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignal[i][j]->SetName(Form("%sOS_DYsignal_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
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

    float Wjet_norm = MC_nEvents_Wjet->Integral();
    massPeakOS_Wjet[i][j]->Scale(s.Wjet_XS/Wjet_norm);

    float ttbar_norm = MC_nEvents_ttbar->Integral();
    massPeakOS_ttbar[i][j]->Scale(s.ttbar_XS/ttbar_norm);

    //remove entries not passing photon selection or SS from signal
    //data
    massPeakOS_minusSSAndPhoton[i][j] = (TH1D*)massPeakOS[i][j]->Clone(Form("%sOS_minusSSAndPhoton_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_minusSSAndPhoton[i][j]->Add( massPeakOS_photons[i][j] , -1);
    massPeakOS_minusSSAndPhoton[i][j]->Add( massPeakSS[i][j], -1);
    
    //same for DY sample (only photon here)
    massPeakOS_DYsignalMinusPhoton[i][j] = (TH1D*)massPeakOS_DYsignal[i][j]->Clone(Form("%sOS_DYsignalMinusPhoton_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhoton[i][j]->Add(massPeakOS_DYphoton[i][j],-1);

    //add the other backgrounds
    massPeakOS_DYsignalMinusPhotonPlusBkg[i][j] = (TH1D*)massPeakOS_DYsignalMinusPhoton[i][j]->Clone(Form("%sOS_DYsignalMinusPhotonPlusBkg_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
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
    

    //normalize MC fractions to data with a ratio of counts between subtracted MC data
    //tau
    bkg_tau[i][j] = (TH1D*) fraction_tau[i][j]->Clone(Form("%sBkg_tau_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_tau[i][j]->Multiply(massPeakOS_minusSSAndPhoton[i][j]);
    
    //ttbar
    bkg_ttbar[i][j] = (TH1D*) fraction_ttbar[i][j]->Clone(Form("%sBkg_ttbar_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_ttbar[i][j]->Multiply(massPeakOS_minusSSAndPhoton[i][j]);

    //Wjet
    bkg_Wjet[i][j] = (TH1D*) fraction_Wjet[i][j]->Clone(Form("%sBkg_Wjet_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_Wjet[i][j]->Multiply(massPeakOS_minusSSAndPhoton[i][j]);

    //remove these backgrounds
    massPeakOS_minusAll[i][j] = (TH1D*) massPeakOS_minusSSAndPhoton[i][j]->Clone(Form("%sOS_minusAll_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_minusAll[i][j]->Add( bkg_tau[i][j], -1);
    massPeakOS_minusAll[i][j]->Add( bkg_ttbar[i][j], -1);
    massPeakOS_minusAll[i][j]->Add( bkg_Wjet[i][j], -1);

    massPeakOS[i][j]->Write();
    massPeakSS[i][j]->Write();
    massPeakOS_photons[i][j]->Write();
    massPeakOS_DYsignalMinusPhoton[i][j]->Write();
    massPeakOS_minusSSAndPhoton[i][j]->Write();
    massPeakOS_minusAll[i][j]->Write();
    fraction_Wjet[i][j]->Write(); 
    fraction_ttbar[i][j]->Write(); 
    fraction_tau[i][j]->Write(); 
    bkg_Wjet[i][j]->Write(); 
    bkg_ttbar[i][j]->Write(); 
    bkg_tau[i][j]->Write(); 
    massPeakOS_DYsignal[i][j]->Write();
    massPeakOS_DYphoton[i][j]->Write();
    massPeakOS_DYtautau[i][j]->Write();
    massPeakOS_Wjet[i][j]->Write();
    massPeakOS_ttbar[i][j]->Write();

    //Draw the histogram Stack
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetBottomMargin(0.2);
    
    if(j==0) massPeakOS[i][j]->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
    if(j==1) massPeakOS[i][j]->GetXaxis()->SetTitle("p_{T} (GeV)");
    if(j==2) massPeakOS[i][j]->GetXaxis()->SetTitle("y");
    if(j==3) massPeakOS[i][j]->GetXaxis()->SetTitle("yield");
    
    if(j==0 || j==3) massPeakOS[i][j]->GetYaxis()->SetTitle("Entries");
    if(j==1) massPeakOS[i][j]->GetYaxis()->SetTitle("#frac{dN_{Z}}{dp_{T}} (GeV^{-1})");
    if(j==2) massPeakOS[i][j]->GetYaxis()->SetTitle("#frac{dN_{Z}}{dy}");
    massPeakOS[i][j]->GetXaxis()->CenterTitle();
    massPeakOS[i][j]->GetYaxis()->CenterTitle();
    massPeakOS[i][j]->SetFillColor(kOrange+1);
    massPeakOS[i][j]->SetMarkerStyle(8);
    massPeakOS[i][j]->SetMarkerColor(kBlack);
    massPeakOS[i][j]->SetLineColor(kBlack);
    massPeakOS[i][j]->SetStats(0);

    TH1D * dummy;
    if(j==1){
      dummy = new TH1D("dummy",";p_{T} (GeV); #frac{dN_{Z}}{dp_{T}} (GeV^{-1})",2,0.1,200);
      dummy->SetBinContent(1,massPeakOS[i][j]->GetMaximum());
      dummy->SetBinContent(2,massPeakOS[i][j]->GetMinimum());
      dummy->SetLineColor(kWhite);
      dummy->SetMarkerColor(kWhite);
      dummy->SetStats(0);
      dummy->Draw();
    }

    massPeakOS[i][j]->Draw("HIST same");

    bkg_ttbar[i][j]->SetFillColor(kGray);
    bkg_ttbar[i][j]->SetLineColor(kBlack);

    bkg_Wjet[i][j]->Add(bkg_ttbar[i][j]);
    bkg_Wjet[i][j]->SetFillColor(kViolet+1);
    bkg_Wjet[i][j]->SetLineColor(kBlack);
    
    bkg_tau[i][j]->Add(bkg_Wjet[i][j]);
    bkg_tau[i][j]->SetFillColor(kGreen+1);
    bkg_tau[i][j]->SetLineColor(kBlack);
    
    massPeakOS_photons[i][j]->Add(bkg_tau[i][j]);
    massPeakOS_photons[i][j]->SetFillColor(kRed+1);
    massPeakOS_photons[i][j]->SetLineColor(kBlack);
    
    massPeakSS[i][j]->Add(massPeakOS_photons[i][j]);
    massPeakSS[i][j]->SetFillColor(kBlue);
    massPeakSS[i][j]->SetLineColor(kBlack); 

    massPeakSS[i][j]->Draw("HIST same");
    massPeakOS_photons[i][j]->Draw("HIST same");
    bkg_tau[i][j]->Draw("HIST same");
    bkg_Wjet[i][j]->Draw("HIST same");
    bkg_ttbar[i][j]->Draw("HIST same");
    massPeakOS[i][j]->Draw("p same");

    TLegend *leg = new TLegend(0.65,0.7,0.9,0.875);
    leg->AddEntry(massPeakOS[i][j],Form("Data (%d-%d%%)",c.getCentBinLow(i),c.getCentBinHigh(i)),"p");
    if(isMu) leg->AddEntry(massPeakOS[i][j],"Z #rightarrow #mu^{+}#mu^{-}","f");
    if(!isMu) leg->AddEntry(massPeakOS[i][j],"Z #rightarrow e^{+}e^{-}","f");
    leg->AddEntry(massPeakSS[i][j],"Same Sign (QCD)","f");
    leg->AddEntry(massPeakOS_photons[i][j],"EM background","f");
    leg->AddEntry(bkg_tau[i][j],"Z #rightarrow #tau^{+}#tau^{-}","f");
    leg->AddEntry(bkg_Wjet[i][j],"W^{#pm} + X","f");
    leg->AddEntry(bkg_ttbar[i][j],"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    if(j==1){
      c1->SetLogx();
      c1->SetLogy();
      dummy->GetXaxis()->SetRangeUser(0.1,199);
      dummy->GetYaxis()->SetRangeUser(0.01,200000);
    }

    c1->RedrawAxis();
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));

    if(j==0) massPeakOS[i][j]->GetYaxis()->SetRangeUser(0.01, massPeakOS[i][j]->GetMaximum()*50);
    if(j==2) massPeakOS[i][j]->GetYaxis()->SetRangeUser(0.3, massPeakOS[i][j]->GetMaximum()*300);
    c1->SetLogy();
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
   
    if(j==1) delete dummy; 
    delete c1;
    delete leg;
  }
  }

  if(isMu)



  out->Close();
  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 7)
  {
    std::cout << "Usage: massPeakPlots_BkgSub <Data file> <DY File> <TTbar File> <WJet File> <isMu> <outTag>" << std::endl;
    return 1;
  }  

  std::string data = argv[1];
  std::string DY = argv[2];
  std::string ttbar = argv[3];
  std::string Wjet = argv[4];
  bool isMu = (bool)std::atoi(argv[5]);
  std::string outTag = argv[6];
   
  plotMassPeaks_BkgSub(data, DY, ttbar, Wjet, isMu, outTag);
  return 0; 
}

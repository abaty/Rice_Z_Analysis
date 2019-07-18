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

void plotMassPeaks_BkgSub(std::string data_, std::string DY_, std::string ttbar_, std::string Wjet_, bool isMu){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TH1D * massPeakOS_photons[nBins]; 
  TH1D * massPeakOS_minusSSAndPhoton[nBins];
  TH1D * massPeakOS_minusAll[nBins];

  TH1D * massPeakOS_DYsignal[nBins]; 
  TH1D * massPeakOS_DYphoton[nBins]; 
  TH1D * massPeakOS_DYtautau[nBins]; 
  TH1D * massPeakOS_DYsignalMinusPhoton[nBins];

  TH1D * massPeakOS_Wjet[nBins]; 
  TH1D * massPeakOS_ttbar[nBins]; 
 
  TH1D * MC_nEvents_DY;
  TH1D * MC_nEvents_Wjet;
  TH1D * MC_nEvents_ttbar;

  TH1D * fraction_tau[nBins];
  TH1D * fraction_Wjet[nBins];
  TH1D * fraction_ttbar[nBins];

  TH1D * bkg_tau[nBins];
  TH1D * bkg_Wjet[nBins];
  TH1D * bkg_ttbar[nBins];
 
  TFile * data = TFile::Open(data_.c_str(),"read");
  TFile * DY = TFile::Open(DY_.c_str(),"read");
  TFile * ttbar = TFile::Open(ttbar_.c_str(),"read");
  TFile * Wjet = TFile::Open(Wjet_.c_str(),"read");
  
  MC_nEvents_DY = (TH1D*)DY->Get("nEvents");  
  MC_nEvents_Wjet = (TH1D*)Wjet->Get("nEvents");  
  MC_nEvents_ttbar= (TH1D*)ttbar->Get("nEvents");  

  TFile * out = TFile::Open(Form("backgroundSubtraction_isMu%d.root",(int)isMu),"recreate");
  
  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = (TH1D*)data->Get(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS[i] = (TH1D*)data->Get(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_photons[i] = (TH1D*)data->Get(Form("massPeakOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS_DYsignal[i] = (TH1D*)DY->Get(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignal[i]->SetName(Form("massPeakOS_DYsignal_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphoton[i] = (TH1D*)DY->Get(Form("massPeakOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphoton[i]->SetName(Form("massPeakOS_DY_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYtautau[i] = (TH1D*)DY->Get(Form("massPeakOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS_Wjet[i] = (TH1D*)Wjet->Get(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbar[i] = (TH1D*)ttbar->Get(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    massPeakOS_Wjet[i]->SetName(Form("massPeakOS_Wjet_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbar[i]->SetName(Form("massPeakOS_ttbar_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
 
    //do normalizations for MC
    float DY_norm = MC_nEvents_DY->Integral();
    massPeakOS_DYsignal[i]->Scale(s.DY_XS/DY_norm);
    massPeakOS_DYphoton[i]->Scale(s.DY_XS/DY_norm);
    massPeakOS_DYtautau[i]->Scale(s.DY_XS/DY_norm);

    float Wjet_norm = MC_nEvents_Wjet->Integral();
    massPeakOS_Wjet[i]->Scale(s.Wjet_XS/Wjet_norm);

    float ttbar_norm = MC_nEvents_ttbar->Integral();
    massPeakOS_ttbar[i]->Scale(s.ttbar_XS/ttbar_norm);

    //remove entries not passing photon selection or SS from signal
    //data
    massPeakOS_minusSSAndPhoton[i] = (TH1D*)massPeakOS[i]->Clone(Form("massPeakOS_minusSSAndPhoton_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_minusSSAndPhoton[i]->Add( massPeakOS_photons[i] , -1);
    massPeakOS_minusSSAndPhoton[i]->Add( massPeakSS[i], -1);
    
    //same for DY sample (only photon here)
    massPeakOS_DYsignalMinusPhoton[i] = (TH1D*)massPeakOS_DYsignal[i]->Clone(Form("massPeakOS_DYsignalMinusPhoton_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhoton[i]->Add(massPeakOS_DYphoton[i],-1);

    //Calculate the fractions for tautau, Wjet, TTbar
    //tautau
    fraction_tau[i] = (TH1D*)massPeakOS_DYtautau[i]->Clone(Form("fraction_tau_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_tau[i]->Divide( massPeakOS_DYsignalMinusPhoton[i] );

    //ttbar
    fraction_ttbar[i] = (TH1D*)massPeakOS_ttbar[i]->Clone(Form("fraction_ttbar_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_ttbar[i]->Divide( massPeakOS_DYsignalMinusPhoton[i] );
    
    //Wjet
    fraction_Wjet[i] = (TH1D*)massPeakOS_Wjet[i]->Clone(Form("fraction_Wjet_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_Wjet[i]->Divide( massPeakOS_DYsignalMinusPhoton[i] );
    

    //normalize MC fractions to data with a ratio of counts between subtracted MC data
    //tau
    bkg_tau[i] = (TH1D*) fraction_tau[i]->Clone(Form("bkg_tau_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_tau[i]->Multiply(massPeakOS_minusSSAndPhoton[i]);
    
    //ttbar
    bkg_ttbar[i] = (TH1D*) fraction_ttbar[i]->Clone(Form("bkg_ttbar_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_ttbar[i]->Multiply(massPeakOS_minusSSAndPhoton[i]);

    //Wjet
    bkg_Wjet[i] = (TH1D*) fraction_Wjet[i]->Clone(Form("bkg_Wjet_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    bkg_Wjet[i]->Multiply(massPeakOS_minusSSAndPhoton[i]);

    //remove these backgrounds
    massPeakOS_minusAll[i] = (TH1D*) massPeakOS_minusSSAndPhoton[i]->Clone(Form("massPeakOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_minusAll[i]->Add( bkg_tau[i], -1);
    massPeakOS_minusAll[i]->Add( bkg_ttbar[i], -1);
    massPeakOS_minusAll[i]->Add( bkg_Wjet[i], -1);

    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    massPeakOS_photons[i]->Write();
    massPeakOS_DYsignalMinusPhoton[i]->Write();
    massPeakOS_minusSSAndPhoton[i]->Write();
    massPeakOS_minusAll[i]->Write();
    fraction_Wjet[i]->Write(); 
    fraction_ttbar[i]->Write(); 
    fraction_tau[i]->Write(); 
    bkg_Wjet[i]->Write(); 
    bkg_ttbar[i]->Write(); 
    bkg_tau[i]->Write(); 
    massPeakOS_DYsignal[i]->Write();
    massPeakOS_DYphoton[i]->Write();
    massPeakOS_DYtautau[i]->Write();
    massPeakOS_Wjet[i]->Write();
    massPeakOS_ttbar[i]->Write();

    //Draw the histogram Stack
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetErrorX(0);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetBottomMargin(0.2);
    
    massPeakOS[i]->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
    massPeakOS[i]->GetYaxis()->SetTitle("Entries");
    massPeakOS[i]->GetXaxis()->CenterTitle();
    massPeakOS[i]->GetYaxis()->CenterTitle();
    massPeakOS[i]->SetFillColor(kOrange+1);
    massPeakOS[i]->SetMarkerStyle(8);
    massPeakOS[i]->SetMarkerColor(kBlack);
    massPeakOS[i]->SetLineColor(kBlack);
    massPeakOS[i]->SetStats(0);
    massPeakOS[i]->Draw("HIST");

    bkg_ttbar[i]->SetFillColor(kGray);
    bkg_ttbar[i]->SetLineColor(kBlack);

    bkg_Wjet[i]->Add(bkg_ttbar[i]);
    bkg_Wjet[i]->SetFillColor(kViolet+1);
    bkg_Wjet[i]->SetLineColor(kBlack);
    
    bkg_tau[i]->Add(bkg_Wjet[i]);
    bkg_tau[i]->SetFillColor(kGreen+1);
    bkg_tau[i]->SetLineColor(kBlack);
    
    massPeakOS_photons[i]->Add(bkg_tau[i]);
    massPeakOS_photons[i]->SetFillColor(kRed+1);
    massPeakOS_photons[i]->SetLineColor(kBlack);
    
    massPeakSS[i]->Add(massPeakOS_photons[i]);
    massPeakSS[i]->SetFillColor(kBlue);
    massPeakSS[i]->SetLineColor(kBlack); 

    massPeakSS[i]->Draw("HIST same");
    massPeakOS_photons[i]->Draw("HIST same");
    bkg_tau[i]->Draw("HIST same");
    bkg_Wjet[i]->Draw("HIST same");
    bkg_ttbar[i]->Draw("HIST same");
    massPeakOS[i]->Draw("p same");

    TLegend *leg = new TLegend(0.625,0.7,0.875,0.875);
    leg->AddEntry(massPeakOS[i],Form("Data (%d-%d%%)",c.getCentBinLow(i),c.getCentBinHigh(i)),"p");
    if(isMu) leg->AddEntry(massPeakOS[i],"Z #rightarrow #mu^{+}#mu^{-}","f");
    if(!isMu) leg->AddEntry(massPeakOS[i],"Z #rightarrow e^{+}e^{-}","f");
    leg->AddEntry(massPeakSS[i],"Same Sign (QCD)","f");
    leg->AddEntry(massPeakOS_photons[i],"EM background","f");
    leg->AddEntry(bkg_tau[i],"Z #rightarrow #tau^{+}#tau^{-}","f");
    leg->AddEntry(bkg_Wjet[i],"W^{#pm} + X","f");
    leg->AddEntry(bkg_ttbar[i],"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    c1->RedrawAxis();
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d.png",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d.pdf",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d.C",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS[i]->GetYaxis()->SetRangeUser(0.01, massPeakOS[i]->GetMaximum()*50);
    c1->SetLogy();
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d_log.png",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d_log.pdf",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks_withBkgSub/massPeak_withSub_isMu%d_%d_%d_log.C",(int)isMu, c.getCentBinLow(i),c.getCentBinHigh(i)));
    


  }

  if(isMu)



  out->Close();
  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 6)
  {
    std::cout << "Usage: massPeakPlots_BkgSub <Data file> <DY File> <TTbar File> <WJet File> <isMu>" << std::endl;
    return 1;
  }  

  std::string data = argv[1];
  std::string DY = argv[2];
  std::string ttbar = argv[3];
  std::string Wjet = argv[4];
  bool isMu = (bool)std::atoi(argv[5]);
   
  plotMassPeaks_BkgSub(data, DY, ttbar, Wjet, isMu);
  return 0; 
}

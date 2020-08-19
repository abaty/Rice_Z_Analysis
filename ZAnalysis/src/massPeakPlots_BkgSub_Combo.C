#include "include/RooUnfold/RooUnfoldResponse.h"
#include "include/RooUnfold/RooUnfoldBayes.h"
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

void plotMassPeaks_BkgSub(std::string data_, std::string DY_, std::string ttbar_, std::string Wjet_, std::string dataE_, std::string DYE_, std::string ttbarE_, std::string WjetE_){
  gErrorIgnoreLevel = kWarning;
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
  
  TH1D * massPeakOSE[nBins][4][5]; 
  TH1D * massPeakSSE[nBins][4][5]; 
  TH1D * massPeakSS_ChargeFlipCorrectedE[nBins][4][5]; 
  TH1D * massPeakOS_photonsE[nBins][4][5]; 
  TH1D * massPeakOS_minusSSAndPhotonE[nBins][4][5];
  TH1D * massPeakOS_minusAllE[nBins][4][5];

  TH1D * massPeakOS_DYsignalE[nBins][4]; 
  TH1D * massPeakOS_DYphotonE[nBins][4]; 
  TH1D * massPeakOS_DYtautauE[nBins][4]; 
  TH1D * massPeakOS_DYsignalMinusPhotonE[nBins][4];
  TH1D * massPeakOS_DYsignalMinusPhotonPlusSSE[nBins][4];
  TH1D * massPeakOS_DYsignalMinusPhotonPlusBkgE[nBins][4];
  TH1D * massPeakSS_DYE[nBins][4];

  TH1D * massPeakOS_WjetE[nBins][4]; 
  TH1D * massPeakOS_ttbarE[nBins][4]; 
 
  TH1D * MC_nEvents_DYE;
  TH1D * MC_nEvents_WjetE;
  TH1D * MC_nEvents_ttbarE;

  TH1D * fraction_tauE[nBins][4];
  TH1D * fraction_WjetE[nBins][4];
  TH1D * fraction_ttbarE[nBins][4];

  TH1D * bkg_tauE[nBins][4][5];
  TH1D * bkg_WjetE[nBins][4][5];
  TH1D * bkg_ttbarE[nBins][4][5];
 
  TFile * data = TFile::Open(data_.c_str(),"read");
  TFile * DY = TFile::Open(DY_.c_str(),"read");
  TFile * ttbar = TFile::Open(ttbar_.c_str(),"read");
  TFile * Wjet = TFile::Open(Wjet_.c_str(),"read");
  
  TFile * dataE = TFile::Open(dataE_.c_str(),"read");
  TFile * DYE = TFile::Open(DYE_.c_str(),"read");
  TFile * ttbarE = TFile::Open(ttbarE_.c_str(),"read");
  TFile * WjetE = TFile::Open(WjetE_.c_str(),"read");
  
  MC_nEvents_DY = (TH1D*)DY->Get("nEvents");  
  MC_nEvents_Wjet = (TH1D*)Wjet->Get("nEvents");  
  MC_nEvents_ttbar= (TH1D*)ttbar->Get("nEvents");  
  
  MC_nEvents_DYE = (TH1D*)DYE->Get("nEvents");  
  MC_nEvents_WjetE = (TH1D*)WjetE->Get("nEvents");  
  MC_nEvents_ttbarE= (TH1D*)ttbarE->Get("nEvents");  

  for(int i = 0; i<nBins; i++){
  if(i!=11) continue;//0-100%
  for(int j = 0; j<4; j++){
    if(j!=0) continue;
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
    massPeakOS_DYsignal[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakOS_DYphoton[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakOS_DYtautau[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakSS_DY[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);

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
      //bkg_tau[i][j][k] = (TH1D*) fraction_tau[i][j]->Clone(Form("%sBkg_tau%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_tau[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);
      bkg_tau[i][j][k] = (TH1D*) massPeakOS_DYtautau[i][j]->Clone(Form("%sBkg_tau%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_tau[i][j][k]->Scale(s.muLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9 * 67.6/70.0); //5.649 is 0-100% TAA, 7644 is the total MB xsection for PbPb in mb, last fraction is because we used older Ncol scaling values that had sigma of 70 for MB pp
    
      //ttbar
      //bkg_ttbar[i][j][k] = (TH1D*) fraction_ttbar[i][j]->Clone(Form("%sBkg_ttbar%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_ttbar[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);
      bkg_ttbar[i][j][k] = (TH1D*) massPeakOS_ttbar[i][j]->Clone(Form("%sBkg_ttbar%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_ttbar[i][j][k]->Scale(s.muLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9 * 67.6/70.0);//see tau tau comment above

      //Wjet
      //bkg_Wjet[i][j][k] = (TH1D*) fraction_Wjet[i][j]->Clone(Form("%sBkg_Wjet%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_Wjet[i][j][k]->Multiply(massPeakOS_minusSSAndPhoton[i][j][k]);
      bkg_Wjet[i][j][k] = (TH1D*) massPeakOS_Wjet[i][j]->Clone(Form("%sBkg_Wjet%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_Wjet[i][j][k]->Scale(s.muLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9* 67.6/70.0);//see W jet comment

    //remove these backgrounds
      massPeakOS_minusAll[i][j][k] = (TH1D*) massPeakOS_minusSSAndPhoton[i][j][k]->Clone(Form("%sOS_minusAll%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_minusAll[i][j][k]->Add( bkg_tau[i][j][k], -1);
      massPeakOS_minusAll[i][j][k]->Add( bkg_ttbar[i][j][k], -1);
      massPeakOS_minusAll[i][j][k]->Add( bkg_Wjet[i][j][k], -1);
    }
  
  //repeat for electrons
  for(int j = 0; j<4; j++){
    if(j!=0) continue;
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      massPeakOSE[i][j][k] = (TH1D*)dataE->Get(Form("%sOS_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakSSE[i][j][k] = (TH1D*)dataE->Get(Form("%sSS_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_photonsE[i][j][k] = (TH1D*)dataE->Get(Form("%sOS_ptLT0p5acoLT0p001_withEff%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    }

    massPeakOS_DYsignalE[i][j] = (TH1D*)DYE->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalE[i][j]->SetName(Form("%sOS_DYsignal_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_DYE[i][j] = (TH1D*)DYE->Get(Form("%sSS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_DYE[i][j]->SetName(Form("%sSS_DYsignal_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphotonE[i][j] = (TH1D*)DYE->Get(Form("%sOS_ptLT0p5acoLT0p001_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYphotonE[i][j]->SetName(Form("%sOS_DY_ptLT0p5acoLT0p001_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYtautauE[i][j] = (TH1D*)DYE->Get(Form("%sOS_TauTau_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));

    massPeakOS_WjetE[i][j] = (TH1D*)WjetE->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbarE[i][j] = (TH1D*)ttbarE->Get(Form("%sOS_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    massPeakOS_WjetE[i][j]->SetName(Form("%sOS_Wjet_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_ttbarE[i][j]->SetName(Form("%sOS_ttbar_withEff_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
 
    //do normalizations for MC
    float DY_norm = MC_nEvents_DYE->Integral();
    massPeakOS_DYsignalE[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakOS_DYphotonE[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakOS_DYtautauE[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);
    massPeakSS_DYE[i][j]->Scale(s.DY_XS*s.crossSectionModifier/DY_norm);

    float Wjet_norm = MC_nEvents_WjetE->Integral();
    massPeakOS_WjetE[i][j]->Scale(s.Wjet_XS/Wjet_norm);

    float ttbar_norm = MC_nEvents_ttbarE->Integral();
    massPeakOS_ttbarE[i][j]->Scale(s.ttbar_XS/ttbar_norm);


    //Remove photon background from signal (these get cut in data as well)
    massPeakOS_DYsignalMinusPhotonE[i][j] = (TH1D*)massPeakOS_DYsignalE[i][j]->Clone(Form("%sOS_DYsignalMinusPhoton_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhotonE[i][j]->Add(massPeakOS_DYphotonE[i][j],-1);
    //add back the same sign pairs that we lose from charge-flipping (note that I think we should not do this now)
    massPeakOS_DYsignalMinusPhotonPlusSSE[i][j] = (TH1D*)massPeakOS_DYsignalMinusPhotonE[i][j]->Clone(Form("%sOS_DYsignalMinusPhotonPlusSS_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    //massPeakOS_DYsignalMinusPhotonPlusSS[i][j]->Add(massPeakSS_DY[i][j],1);

    //add the other backgrounds
    massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j] = (TH1D*)massPeakOS_DYsignalMinusPhotonPlusSSE[i][j]->Clone(Form("%sOS_DYsignalMinusPhotonPlusBkg_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j]->Add(massPeakOS_DYtautauE[i][j]);
    massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j]->Add(massPeakOS_ttbarE[i][j]);
    massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j]->Add(massPeakOS_WjetE[i][j]);

    //Calculate the fractions for tautau, Wjet, TTbar
    //tautau
    fraction_tauE[i][j] = (TH1D*)massPeakOS_DYtautauE[i][j]->Clone(Form("%sFraction_tau_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_tauE[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j] );

    //ttbar
    fraction_ttbarE[i][j] = (TH1D*)massPeakOS_ttbarE[i][j]->Clone(Form("%sFraction_ttbar_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_ttbarE[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j] );
    
    //Wjet
    fraction_WjetE[i][j] = (TH1D*)massPeakOS_WjetE[i][j]->Clone(Form("%sFraction_Wjet_%d_%d",h.name.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
    fraction_WjetE[i][j]->Divide( massPeakOS_DYsignalMinusPhotonPlusBkgE[i][j] );
   
    //remove entries not passing photon selection or SS from signal
    //data
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      massPeakOS_minusSSAndPhotonE[i][j][k] = (TH1D*)massPeakOSE[i][j][k]->Clone(Form("%sOS_minusSSAndPhoton%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_minusSSAndPhotonE[i][j][k]->Add( massPeakOS_photonsE[i][j][k] , -1);
      
      //take charge flip into account
      massPeakSS_ChargeFlipCorrectedE[i][j][k] = (TH1D*) massPeakSSE[i][j][k]->Clone(Form("%sSS_ChargeFlipCorrected%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i))); 
      //scan along bins looking for nonzero counts, and adjust based on MC estimation
      //commented lines are for adding the charge swaps back as signal (our MC efficiency assumes OS reco, so the same signs should not be added back in I guess...)
      for(int b = 0; b<massPeakSS_ChargeFlipCorrectedE[i][j][k]->GetSize(); b++){
        float SSYield = massPeakSS_ChargeFlipCorrectedE[i][j][k]->GetBinContent(b);
        if( SSYield == 0 ) continue;
        
        float estimate = ( massPeakSS_DYE[i][j]->GetBinContent(b) / massPeakOS_DYsignalMinusPhotonPlusSSE[i][j]->GetBinContent(b) ) * massPeakOS_minusSSAndPhotonE[i][j][k]->GetBinContent(b);
        if( SSYield > estimate ){
          massPeakSS_ChargeFlipCorrectedE[i][j][k]->SetBinContent(b, SSYield-estimate);
         // massPeakOS_minusSSAndPhoton[i][j][k]->SetBinContent(b, massPeakOS_minusSSAndPhoton[i][j][k]->GetBinContent(b) + estimate );
        }
        else{
          massPeakSS_ChargeFlipCorrectedE[i][j][k]->SetBinContent(b, 0);
         // massPeakOS_minusSSAndPhoton[i][j][k]->SetBinContent(b, massPeakOS_minusSSAndPhoton[i][j][k]->GetBinContent(b) + SSYield );
        }
      }   
 
      //remove remaining SS background
      massPeakOS_minusSSAndPhotonE[i][j][k]->Add( massPeakSS_ChargeFlipCorrectedE[i][j][k], -1);
    }    

    //normalize MC fractions to data with a ratio of counts between subtracted MC data
    //tau
    for(int k = 0; k<5; k++){
      if(j==0 && k>0) continue;
      //bkg_tauE[i][j][k] = (TH1D*) fraction_tauE[i][j]->Clone(Form("%sBkg_tau%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_tauE[i][j][k]->Multiply(massPeakOS_minusSSAndPhotonE[i][j][k]);
      bkg_tauE[i][j][k] = (TH1D*) massPeakOS_DYtautauE[i][j]->Clone(Form("%sBkg_tau%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_tauE[i][j][k]->Scale(s.eLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9* 67.6/70.0); //5.649 is 0-100% TAA, 7644 is the total MB xsection for PbPb in mb, last fraction is because we used older Ncol scaling values that had sigma of 70 for MB pp
    
      //ttbar
      //bkg_ttbarE[i][j][k] = (TH1D*) fraction_ttbarE[i][j]->Clone(Form("%sBkg_ttbar%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_ttbarE[i][j][k]->Multiply(massPeakOS_minusSSAndPhotonE[i][j][k]);
      bkg_ttbarE[i][j][k] = (TH1D*) massPeakOS_ttbarE[i][j]->Clone(Form("%sBkg_ttbar%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_ttbarE[i][j][k]->Scale(s.eLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9 * 67.6/70.0);//see tau tau comment above

      //Wjet
      //bkg_WjetE[i][j][k] = (TH1D*) fraction_WjetE[i][j]->Clone(Form("%sBkg_Wjet%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      //bkg_WjetE[i][j][k]->Multiply(massPeakOS_minusSSAndPhotonE[i][j][k]);
      bkg_WjetE[i][j][k] = (TH1D*) massPeakOS_WjetE[i][j]->Clone(Form("%sBkg_Wjet%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      bkg_WjetE[i][j][k]->Scale(s.eLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9* 67.6/70.0);//see W jet comment

    //remove these backgrounds
      massPeakOS_minusAllE[i][j][k] = (TH1D*) massPeakOS_minusSSAndPhotonE[i][j][k]->Clone(Form("%sOS_minusAll%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));
      massPeakOS_minusAllE[i][j][k]->Add( bkg_tauE[i][j][k], -1);
      massPeakOS_minusAllE[i][j][k]->Add( bkg_ttbarE[i][j][k], -1);
      massPeakOS_minusAllE[i][j][k]->Add( bkg_WjetE[i][j][k], -1);
    }


    //Draw the histogram Stack
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    //gStyle->SetErrorX(0);
    TCanvas * c1 = new TCanvas("c1","c1",1600,800);
    c1->SetBorderSize(0);
    TPad * p1 = new TPad("p1","p1",0,0.24,0.53,1,0);
    TPad * p2 = new TPad("p2","p2",0,0,0.53,0.24,0);
    TPad * p3 = new TPad("p3","p3",0.53,0.24,1,1,0);
    TPad * p4 = new TPad("p4","p4",0.53,0,1,0.24,0);
    c1->SetLineWidth(0);
    p1->SetBottomMargin(0);
    p1->SetLeftMargin(0.15);
    p1->SetRightMargin(0);
    p1->SetTopMargin(0.1);
    p1->SetBorderSize(0);
    p1->Draw();
    p2->SetTopMargin(0);
    p2->SetLeftMargin(0.15);
    p2->SetRightMargin(0.0);
    p2->SetBottomMargin(0.5);
    p2->SetBorderSize(0);
    p2->Draw();
    p3->SetBottomMargin(0);
    p3->SetLeftMargin(0);
    p3->SetRightMargin(0.04);
    p3->SetTopMargin(0.1);
    p3->SetBorderSize(0);
    p3->Draw();
    p4->SetTopMargin(0);
    p4->SetLeftMargin(0);
    p4->SetRightMargin(0.04);
    p4->SetBottomMargin(0.5);
    p4->SetBorderSize(0);
    p4->Draw();
    p1->cd();
    
    massPeakOS[i][j][0]->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
    massPeakOSE[i][j][0]->GetXaxis()->SetTitle("m_{ee} (GeV)");
    
    if(j==0 || j==3) massPeakOS[i][j][0]->GetYaxis()->SetTitle("Entries");
    massPeakOS[i][j][0]->GetXaxis()->CenterTitle();
    massPeakOSE[i][j][0]->GetXaxis()->CenterTitle();
    massPeakOS[i][j][0]->GetYaxis()->CenterTitle();
    massPeakOS[i][j][0]->SetFillColor(kOrange+1);
    massPeakOS[i][j][0]->SetMarkerStyle(8);
    massPeakOS[i][j][0]->SetMarkerColor(kBlack);
    massPeakOS[i][j][0]->SetLineColor(kBlack);
    massPeakOS[i][j][0]->SetStats(0);
    
    massPeakOSE[i][j][0]->SetFillColor(kOrange+1);
    massPeakOSE[i][j][0]->SetMarkerStyle(8);
    massPeakOSE[i][j][0]->SetMarkerColor(kBlack);
    massPeakOSE[i][j][0]->SetLineColor(kBlack);
    massPeakOSE[i][j][0]->SetStats(0);

    if(j==0){
      massPeakOS[i][j][0]->GetYaxis()->SetLabelSize(0.07);
      massPeakOS[i][j][0]->GetYaxis()->SetTitleSize(0.09);
      massPeakOS[i][j][0]->GetYaxis()->SetTitleOffset(0.79);
    }
   
    float maxmassOS = massPeakOS[i][j][0]->GetMaximum(); 
    if(j==0) massPeakOS[i][j][0]->GetYaxis()->SetRangeUser(0.15, maxmassOS*5);
    if(j==0) massPeakOSE[i][j][0]->GetYaxis()->SetRangeUser(0.15, maxmassOS*5);

    massPeakOS[i][j][0]->Draw("p same");
    p3->cd();
    massPeakOSE[i][j][0]->Draw("p same");
    p1->cd();

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
    
    bkg_ttbarE[i][j][0]->SetFillColor(kGray);
    bkg_ttbarE[i][j][0]->SetLineColor(kBlack);

    bkg_WjetE[i][j][0]->Add(bkg_ttbarE[i][j][0]);
    bkg_WjetE[i][j][0]->SetFillColor(kViolet+1);
    bkg_WjetE[i][j][0]->SetLineColor(kBlack);
    
    bkg_tauE[i][j][0]->Add(bkg_WjetE[i][j][0]);
    bkg_tauE[i][j][0]->SetFillColor(kGreen+1);
    bkg_tauE[i][j][0]->SetLineColor(kBlack);
    
    massPeakOS_photonsE[i][j][0]->Add(bkg_tauE[i][j][0]);
    massPeakOS_photonsE[i][j][0]->SetFillColor(kRed+1);
    massPeakOS_photonsE[i][j][0]->SetLineColor(kBlack);
    
    massPeakSS_ChargeFlipCorrectedE[i][j][0]->Add(massPeakOS_photonsE[i][j][0]);
    massPeakSS_ChargeFlipCorrectedE[i][j][0]->SetFillColor(kBlue);
    massPeakSS_ChargeFlipCorrectedE[i][j][0]->SetLineColor(kBlack); 

    //normalize our signal to the number of data events minus the ones already accounted for
    //float integralDYSignal = massPeakOS_DYsignalMinusPhoton[i][j]->Integral();
    //float dataIntegral = massPeakOS[i][j][0]->Integral() - massPeakSS_ChargeFlipCorrected[i][j][0]->Integral();
    //massPeakOS_DYsignalMinusPhoton[i][j]->Scale(dataIntegral/integralDYSignal);
    
    massPeakOS_DYsignalMinusPhoton[i][j]->Scale(s.muLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9 * 67.6/70.0);//see tau tau background comment for magic numbers;
    massPeakOS_DYsignalMinusPhoton[i][j]->Add(massPeakSS_ChargeFlipCorrected[i][j][0]);
    massPeakOS_DYsignalMinusPhoton[i][j]->SetFillColor(kOrange+1);
    massPeakOS_DYsignalMinusPhoton[i][j]->SetLineColor(kBlack);
    
    //float integralDYSignalE = massPeakOS_DYsignalMinusPhotonE[i][j]->Integral();
    //float dataIntegralE = massPeakOSE[i][j][0]->Integral() - massPeakSS_ChargeFlipCorrectedE[i][j][0]->Integral();
    //massPeakOS_DYsignalMinusPhotonE[i][j]->Scale(dataIntegralE/integralDYSignalE);
    massPeakOS_DYsignalMinusPhotonE[i][j]->Scale(s.eLumi / s.netLumi * s.Nmb / 0.9 * 5.649e-9* 67.6/70.0);//see tau tau background comment for magic numbers;
    massPeakOS_DYsignalMinusPhotonE[i][j]->Add(massPeakSS_ChargeFlipCorrectedE[i][j][0]);
    massPeakOS_DYsignalMinusPhotonE[i][j]->SetFillColor(kOrange+1);
    massPeakOS_DYsignalMinusPhotonE[i][j]->SetLineColor(kBlack);

    massPeakOS_DYsignalMinusPhoton[i][j]->Draw("HIST same");
    massPeakSS_ChargeFlipCorrected[i][j][0]->Draw("HIST same");
    massPeakOS_photons[i][j][0]->Draw("HIST same");
    bkg_tau[i][j][0]->Draw("HIST same");
    bkg_Wjet[i][j][0]->Draw("HIST same");
    bkg_ttbar[i][j][0]->Draw("HIST same");
    massPeakOS[i][j][0]->Draw("p same");

    p3->cd();
    TLegend *leg = new TLegend(0.65,0.625,0.95,0.85);
    leg->AddEntry(bkg_tau[i][j][0],"Z/#gamma*#rightarrow#tau^{+}#tau^{-}","f");
    leg->AddEntry(bkg_Wjet[i][j][0],"W^{#pm} + X","f");
    leg->AddEntry(bkg_ttbar[i][j][0],"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.065);
    leg->Draw("same");
    
    TLegend *legE = new TLegend(0.03,0.7,0.33,0.85);
    //if(!isMu) leg->AddEntry(massPeakOS_DYsignalMinusPhoton[i][j],"Z #rightarrow e^{+}e^{-}","f");
    legE->AddEntry(massPeakSS_ChargeFlipCorrected[i][j][0],"Same sign (QCD)","f");
    legE->AddEntry(massPeakOS_photons[i][j][0],"EM background","f");
    legE->SetBorderSize(0);
    legE->SetFillStyle(0);
    legE->SetTextSize(0.065);
    legE->Draw("same");
   
    p1->cd();
    TLegend * leg2 = new TLegend(0.65,0.7,0.9,0.85);
    leg2->AddEntry(massPeakOS[i][j][0],"Data","p");
    leg2->AddEntry(massPeakOS_DYsignalMinusPhoton[i][j],"Z/#gamma* MC","f");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.065);
    leg2->Draw("same");
    p3->cd(); 

    massPeakOS_DYsignalMinusPhotonE[i][j]->Draw("HIST same");
    massPeakSS_ChargeFlipCorrectedE[i][j][0]->Draw("HIST same");
    massPeakOS_photonsE[i][j][0]->Draw("HIST same");
    bkg_tauE[i][j][0]->Draw("HIST same");
    bkg_WjetE[i][j][0]->Draw("HIST same");
    bkg_ttbarE[i][j][0]->Draw("HIST same");
    massPeakOSE[i][j][0]->Draw("p same");

    //ratio plot
    TH1D * tempRatio = (TH1D*)massPeakOS[i][j][0]->Clone("tempRatio");
    tempRatio->Divide(massPeakOS_DYsignalMinusPhoton[i][j]);
    tempRatio->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
    p2->cd();
    tempRatio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    tempRatio->GetYaxis()->SetTitleOffset(0.4);
    tempRatio->GetYaxis()->SetRangeUser(0.25,1.75);
    tempRatio->GetYaxis()->SetTitleSize(0.23);
    tempRatio->GetYaxis()->SetTitleOffset(0.29);
    tempRatio->GetYaxis()->SetLabelSize(0.21);
    tempRatio->GetXaxis()->SetTitleSize(0.27);
    tempRatio->GetXaxis()->SetTitleOffset(0.82);
    tempRatio->GetXaxis()->SetLabelSize(0.23);
    tempRatio->GetXaxis()->SetRangeUser(60,120);
    tempRatio->Draw("same");
   
    TBox box = TBox(0,0,0,0);
    box.SetFillColor(kWhite);
    box.SetLineColor(kWhite);
    box.DrawBox(117,-0.35,119.999,0.23);
    TLatex lat = TLatex();
    lat.SetTextSize(0.24);
    lat.SetTextFont(42);
    lat.DrawLatexNDC(.971,0.305,"60");    

    TH1D * tempRatioE = (TH1D*)massPeakOSE[i][j][0]->Clone("tempRatio");
    tempRatioE->Divide(massPeakOS_DYsignalMinusPhotonE[i][j]);
    tempRatioE->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
    p4->cd();
    tempRatioE->GetYaxis()->SetTitle("#frac{Data}{MC}");
    tempRatioE->GetYaxis()->SetTitleOffset(0.4);
    tempRatioE->GetYaxis()->SetRangeUser(0.25,1.75);
    tempRatioE->GetYaxis()->SetTitleSize(0.2);
    tempRatioE->GetYaxis()->SetTitleOffset(0.32);
    tempRatioE->GetYaxis()->SetLabelSize(0.18);
    tempRatioE->GetXaxis()->SetTitleSize(0.27);
    tempRatioE->GetXaxis()->SetTitleOffset(0.82);
    tempRatioE->GetXaxis()->SetLabelSize(0.23);
    tempRatioE->GetXaxis()->SetRangeUser(60,120);
    tempRatioE->Draw("same");
 
    TLine * line1;
    line1 = new TLine(tempRatio->GetXaxis()->GetBinLowEdge(1),1,tempRatio->GetXaxis()->GetBinUpEdge(tempRatio->GetSize()-2),1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");
    p2->cd();
    line1->Draw("same");

    p1->cd();


    p1->RedrawAxis();
    p3->RedrawAxis();

    //c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_Combined.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));
    //c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_Combined.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));
    //c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_Combined.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));

    p1->SetLogy();
    p3->SetLogy();
    CMS_lumi(p3,0,10,1.9,false,true,false, true, true);
    CMS_lumi(p1,0,10,1.9,true,false,true, false, true);

    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log_Combined.png",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log_Combined.pdf",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/%ss_withBkgSub/%s_withSub_isMu%d_%d_%d_log_Combined.C",h.name.at(j).c_str(),h.name.at(j).c_str(),(int)1, c.getCentBinLow(i),c.getCentBinHigh(i)));
   
    delete line1; 
    delete tempRatio;
    delete p1;
    delete p2;
    delete c1;
    delete leg;
  }
  }
  }

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 9)
  {
    std::cout << "Usage: massPeakPlots_BkgSub <Data filemu> <DY Filemu> <TTbar Filemu> <WJet Filemu><Data file e> <DY File e> <TTbar File e> <WJet File e> " << std::endl;
    return 1;
  }  

  std::string data = argv[1];
  std::string DY = argv[2];
  std::string ttbar = argv[3];
  std::string Wjet = argv[4];
  std::string dataE = argv[5];
  std::string DYE = argv[6];
  std::string ttbarE = argv[7];
  std::string WjetE = argv[8];

  plotMassPeaks_BkgSub(data, DY, ttbar, Wjet,dataE, DYE, ttbarE, WjetE);
  return 0; 
}

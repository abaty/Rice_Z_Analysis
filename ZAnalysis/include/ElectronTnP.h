#ifndef ELECTRONTNP
#define ELECTRONTNP

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "TMath.h"

class ElectronTnP{

public:

  ElectronTnP();
  ~ElectronTnP();
  float getZSF(int hiBin, float pt1, float eta, float pt2, float eta2, int idx);

private:

 float getSFFromGraph(float x, TGraphAsymmErrors* g, int idx);
 float getSingleSF(int hiBin, float x, float eta, int idx = 0);
 float getSingleSF_Reco(int hiBin, float x, float eta, int idx = 0);
 float getSingleSF_HLT(float x, float eta,bool isData, int idx = 0);

 TFile * EE_030_f, *EE_30100_f, *EB_030_f, *EB_30100_f;
 TFile * EE_030_fR, *EE_30100_fR, *EB_030_fR, *EB_30100_fR;

 TFile * EE_HLT_Data_f, *EB_HLT_Data_f;
 TFile * EE_HLT_MC_f, *EB_HLT_MC_f;

 TGraphAsymmErrors * EE_HLT_Data;
 TGraphAsymmErrors * EE_HLT_MC;
 TGraphAsymmErrors * EB_HLT_Data;
 TGraphAsymmErrors * EB_HLT_MC;

 TGraphAsymmErrors * EE_030; 
 TGraphAsymmErrors * EE_30100; 
 TGraphAsymmErrors * EB_030; 
 TGraphAsymmErrors * EB_30100;
 
 TGraphAsymmErrors * EE_030R; 
 TGraphAsymmErrors * EE_30100R; 
 TGraphAsymmErrors * EB_030R; 
 TGraphAsymmErrors * EB_30100R;
  
 inline float quad(float x, float y); 
 inline float quad2(float x, float y); 

};

inline float ElectronTnP::quad(float x, float y){
  return TMath::Sqrt(x*x + y*y);
}

inline float ElectronTnP::quad2(float x, float y){
  return x*x + y*y;
}

float ElectronTnP::getSFFromGraph(float x, TGraphAsymmErrors * g, int idx){
  double pt = 0;
  double sf = 1.0;

  for(int i = 0; i<g->GetN(); i++){
    g->GetPoint(i, pt, sf);
  
    if( x >= (pt - g->GetErrorXlow(i)) && x < (pt + g->GetErrorXhigh(i)) ){
      if(idx == 1)  return sf + g->GetErrorYhigh(i);
      if(idx == -1) return sf + g->GetErrorYlow(i);
      return sf;
    }
  }

  return 1.0;
}

float ElectronTnP::getSingleSF(int hiBin, float x, float eta, int idx){
  
  //EB
  if( TMath::Abs(eta)<1.442 ){
    if(hiBin<60){
      return getSFFromGraph(x, EB_030, idx);
    }
    else{
      return getSFFromGraph(x, EB_30100, idx);
    }
  }
  //EE
  else{
    if(hiBin<60){
      return getSFFromGraph(x, EE_030,idx);
    }
    else{
      return getSFFromGraph(x, EE_30100, idx);
    }
  }
  return 1.0;
}

float ElectronTnP::getSingleSF_Reco(int hiBin, float x, float eta, int idx){
  
  //EB
  if( TMath::Abs(eta)<1.442 ){
    if(hiBin<60){
      return getSFFromGraph(x, EB_030R, idx);
    }
    else{
      return getSFFromGraph(x, EB_30100R, idx);
    }
  }
  //EE
  else{
    if(hiBin<60){
      return getSFFromGraph(x, EE_030R,idx);
    }
    else{
      return getSFFromGraph(x, EE_30100R, idx);
    }
  }
  return 1.0;
}

float ElectronTnP::getSingleSF_HLT(float x, float eta, bool isData, int idx){
  
  //EB
  if( TMath::Abs(eta)<1.442 ){
      if(isData) return getSFFromGraph(x, EB_HLT_Data, idx);
      else       return getSFFromGraph(x, EB_HLT_MC, idx);
  }
  //EE
  else{
      if(isData) return getSFFromGraph(x, EE_HLT_Data, idx);
      else       return getSFFromGraph(x, EE_HLT_MC, idx);
  }
  return 1.0;
}

float ElectronTnP::getZSF(int hiBin, float pt1, float eta1, float pt2, float eta2, int idx){

  //for sysetmatic variations eventually;
  idx = idx;

  //if we are >150 GeV, treat it like it's 150 GeV
  if(pt1 >= 150){
    pt1 = 149;
  }

  if(pt2 >= 150){
    pt2 = 149;
  }
  
  //Muon reco scale factors:
  float recoSF1 = getSingleSF_Reco(hiBin, pt1, eta1);
  float recoSF1_var = recoSF1;
  if(idx != 0){
    recoSF1_var = (getSingleSF_Reco(hiBin, pt1, eta1, idx) - recoSF1)/recoSF1;
  }

  float recoSF2 = getSingleSF_Reco(hiBin, pt2, eta2);
  float recoSF2_var = recoSF2;
  if(idx != 0){
    recoSF2_var = (getSingleSF_Reco(hiBin, pt2, eta2, idx) - recoSF2)/recoSF2;
  }

  float recoSF = recoSF1 * recoSF2;
  //if(idx !=0){
  //  if(idx==1 ) recoSF = recoSF * ( 1 + quad(recoSF1_var, recoSF2_var ) );
  //  if(idx==-1) recoSF = recoSF * ( 1 - quad(recoSF1_var, recoSF2_var ) );
  //}

  //Muon ID scale factors:
  float idSF1 = getSingleSF(hiBin, pt1, eta1);
  float idSF1_var = idSF1;
  if(idx != 0){
    idSF1_var = (getSingleSF(hiBin, pt1, eta1, idx) - idSF1)/idSF1;
  }

  float idSF2 = getSingleSF(hiBin, pt2, eta2);
  float idSF2_var = idSF2;
  if(idx != 0){
    idSF2_var = (getSingleSF(hiBin, pt2, eta2, idx) - idSF2)/idSF2;
  }
  

  float idSF = idSF1 * idSF2;
  //if(idx !=0){
  //  if(idx==1 ) idSF = idSF * ( 1 + quad(idSF1_var, idSF2_var ) );
  //  if(idx==-1) idSF = idSF * ( 1 - quad(idSF1_var, idSF2_var ) );
  //}


  //float trigSF1 = 1.0;
  //float trigSF2 = 1.0;

  //float trigSF = trigSF1 * trigSF2;
  //Muon trigger scale factors:
  //efficiency for data (product of inefficiencies because single muon trigger
  float eff1 = getSingleSF_HLT(pt1, eta1, true);
  float eff2 = getSingleSF_HLT(pt2, eta2, true);
  float eff_data = 1 - ((1 - eff1)*(1 - eff2));
  float eff_MC = 1 - ( (1 - getSingleSF_HLT(pt1, eta1, false) ) * (1 - getSingleSF_HLT(pt2, eta2, false)) );

  float trigSF = eff_data/eff_MC;
  float trigSF_varied = trigSF;

  //upward
  if(idx == 1){
    float eff_data_StatU1 = ((1 - ( (1 - getSingleSF_HLT(pt1, eta1, true, 1)) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatU2 = ((1 - ( (1 - eff1) * (1 - getSingleSF_HLT(pt2, eta2, true,1 )) ))- eff_data)/eff_data; 
    float StatU = quad( eff_data_StatU1, eff_data_StatU2); 

    trigSF_varied = StatU;
  }
  //downward
  if(idx == -1){
    float eff_data_StatD1 = ((1 - ( (1 - getSingleSF_HLT(pt1, eta1, true, -1)) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatD2 = ((1 - ( (1 - eff1) * (1 - getSingleSF_HLT(pt2, eta2, true,-1 )) ))- eff_data)/eff_data; 
    float StatD = quad( eff_data_StatD1, eff_data_StatD2); 
    
    trigSF_varied = StatD;
  }

  float SF = recoSF * idSF * trigSF;
  if(idx !=0){
    if(idx==1 ) SF = SF * ( 1 + TMath::Sqrt( quad2(idSF1_var, idSF2_var ) + quad2(recoSF1_var, recoSF2_var) + trigSF_varied*trigSF_varied) );
    if(idx==-1) SF = SF * ( 1 - TMath::Sqrt( quad2(idSF1_var, idSF2_var ) + quad2(recoSF1_var, recoSF2_var) + trigSF_varied*trigSF_varied) );
  }
  return SF;
}

ElectronTnP::ElectronTnP(){
  EE_030_f = TFile::Open("resources/electronSFs/v2_July9/ScaleFactors_PbPb_LooseWP_EE_Centr_0_30_AlpaFixedDataOnly_BWResCBErfExp_preliminary_v2.root","open"); 
  EE_030 = (TGraphAsymmErrors*) EE_030_f->Get("g_scalefactors");

  EE_30100_f = TFile::Open("resources/electronSFs/v2_July9/ScaleFactors_PbPb_LooseWP_EE_Centr_30_100_AlpaFixedDataOnly_BWResCBErfExp_preliminary_v2.root","open"); 
  EE_30100 = (TGraphAsymmErrors*) EE_30100_f->Get("g_scalefactors");
  
  EB_030_f = TFile::Open("resources/electronSFs/v2_July9/ScaleFactors_PbPb_LooseWP_EB_Centr_0_30_AlpaFixedDataOnly_BWResCBErfExp_preliminary_v2.root","open"); 
  EB_030 = (TGraphAsymmErrors*) EB_030_f->Get("g_scalefactors");

  EB_30100_f = TFile::Open("resources/electronSFs/v2_July9/ScaleFactors_PbPb_LooseWP_EB_Centr_30_100_AlpaFixedDataOnly_BWResCBErfExp_preliminary_v2.root","open"); 
  EB_30100 = (TGraphAsymmErrors*) EB_30100_f->Get("g_scalefactors");
 
  //reconstruction 
  EE_030_fR = TFile::Open("resources/electronSFs/RecoSFs_Aug8/Reco_ScaleFactors_PbPb_RECO_EE_Centr_0_30.root","open"); 
  EE_030R = (TGraphAsymmErrors*) EE_030_fR->Get("g_scalefactors");

  EE_30100_fR = TFile::Open("resources/electronSFs/RecoSFs_Aug8/Reco_ScaleFactors_PbPb_RECO_EE_Centr_30_100.root","open"); 
  EE_30100R = (TGraphAsymmErrors*) EE_30100_fR->Get("g_scalefactors");
  
  EB_030_fR = TFile::Open("resources/electronSFs/RecoSFs_Aug8/Reco_ScaleFactors_PbPb_RECO_EB_Centr_0_30.root","open"); 
  EB_030R = (TGraphAsymmErrors*) EB_030_fR->Get("g_scalefactors");

  EB_30100_fR = TFile::Open("resources/electronSFs/RecoSFs_Aug8/Reco_ScaleFactors_PbPb_RECO_EB_Centr_30_100.root","open"); 
  EB_30100R = (TGraphAsymmErrors*) EB_30100_fR->Get("g_scalefactors");

  //trigger
  EE_HLT_Data_f = TFile::Open("resources/electronSFs/HLT_SFs_Aug9/eleTreeEff0_PbPb_LooseWP_EE_Centr_0_100_HLTOnly_Data.root","open");
  EE_HLT_Data = (TGraphAsymmErrors*) EE_HLT_Data_f->Get("Graph");
 
  EE_HLT_MC_f = TFile::Open("resources/electronSFs/HLT_SFs_Aug9/eleTreeEff0_PbPb_LooseWP_EE_Centr_0_100_HLTOnly_MC.root","open");
  EE_HLT_MC = (TGraphAsymmErrors*) EE_HLT_MC_f->Get("Graph");
  
  EB_HLT_Data_f = TFile::Open("resources/electronSFs/HLT_SFs_Aug9/eleTreeEff0_PbPb_LooseWP_EB_Centr_0_100_HLTOnly_Data.root","open");
  EB_HLT_Data = (TGraphAsymmErrors*) EB_HLT_Data_f->Get("Graph");
  
  EB_HLT_MC_f = TFile::Open("resources/electronSFs/HLT_SFs_Aug9/eleTreeEff0_PbPb_LooseWP_EB_Centr_0_100_HLTOnly_MC.root","open");
  EB_HLT_MC = (TGraphAsymmErrors*) EB_HLT_MC_f->Get("Graph");
}

ElectronTnP::~ElectronTnP(){
  EE_030_f->Close();
  EE_30100_f->Close();
  EB_030_f->Close();
  EB_30100_f->Close();
  EE_030_fR->Close();
  EE_30100_fR->Close();
  EB_030_fR->Close();
  EB_30100_fR->Close();

  EE_HLT_Data_f->Close();
  EE_HLT_MC_f->Close();
  EB_HLT_Data_f->Close();
  EB_HLT_MC_f->Close();
}

#endif

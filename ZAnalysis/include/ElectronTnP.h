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

 TFile * EE_030_f, *EE_30100_f, *EB_030_f, *EB_30100_f;

 TGraphAsymmErrors * EE_030; 
 TGraphAsymmErrors * EE_30100; 
 TGraphAsymmErrors * EB_030; 
 TGraphAsymmErrors * EB_30100;
  
 inline float quad(float x, float y); 

};

inline float ElectronTnP::quad(float x, float y){
  return TMath::Sqrt(x*x + y*y);
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
  float recoSF1 = 1.0;
  float recoSF2 = 1.0;

  float recoSF = recoSF1 * recoSF2;


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
  if(idx !=0){
    if(idx==1 ) idSF = idSF * ( 1 + quad(idSF1_var, idSF2_var ) );
    if(idx==-1) idSF = idSF * ( 1 - quad(idSF1_var, idSF2_var ) );
  }


  //Muon trigger scale factors:
  //efficiency for data (product of inefficiencies because single muon trigger
  float trigSF1 = 1.0;
  float trigSF2 = 1.0;

  float trigSF = trigSF1 * trigSF2;

  float SF = recoSF * idSF * trigSF;
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
}

ElectronTnP::~ElectronTnP(){
  EE_030_f->Close();
  EE_30100_f->Close();
  EB_030_f->Close();
  EB_30100_f->Close();
}

#endif

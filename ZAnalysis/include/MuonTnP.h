#ifndef MUONTNP
#define MUONTNP

#include "tnp_weight.h"
#include "TMath.h"

class MuonTnP{

public:
 float getZSF(float pt1, float eta1, float pt2, float eta2, int idx);

private:
 inline float quad2(float x, float y);
 inline float quad2(float x, float y, float z);
 inline float quad(float x, float y);

};

inline float MuonTnP::quad2(float x, float y, float z){
  return x*x+y*y+z*z;
}

inline float MuonTnP::quad2(float x, float y){
  return x*x+y*y;
}

inline float MuonTnP::quad(float x, float y){
  return TMath::Sqrt(quad2(x,y));
}

float MuonTnP::getZSF(float pt1, float eta1, float pt2, float eta2, int idx){
  
  //Muon reco scale factors:
  float recoSF1 = 1.0;
  float recoSF2 = 1.0;

  float recoSF = recoSF1 * recoSF2;

  //******************************************************************************************
  //Muon ID scale factors:
  float idSF1 = tnp_weight_muid_pbpb(eta1 , 0);
  float idSF2 = tnp_weight_muid_pbpb(eta2 , 0);

  float idSF = idSF1 * idSF2;
  float idSF_varied = idSF;
 
  //upward variation
  if(idx == 1){
    //first muon
    float idSF1StatU = (tnp_weight_muid_pbpb(eta1 , 1) - idSF1) / idSF1;// relative stat up uncertainty
    float idSF1SystU = (tnp_weight_muid_pbpb(eta1 , -1) - idSF1) / idSF1;// relative syst up uncertainty
    float idSF1U_2 = quad2(idSF1StatU, idSF1SystU); 
    //second muon
    float idSF2StatU = (tnp_weight_muid_pbpb(eta2 , 1) - idSF2) / idSF2;// relative stat up uncertainty
    float idSF2SystU = (tnp_weight_muid_pbpb(eta2 , -1) - idSF2) / idSF2;// relative syst up uncertainty
    float idSF2U_2 = quad2(idSF2StatU, idSF2SystU);
   
    idSF_varied = idSF1U_2 + idSF2U_2;//relative uncertainty squared
    std::cout << idSF_varied << std::endl;
  }
  //downward variation
  if(idx == -1){
    //first muon
    float idSF1StatD = (tnp_weight_muid_pbpb(eta1 , 2) - idSF1) / idSF1;// relative stat down uncertainty
    float idSF1SystD = (tnp_weight_muid_pbpb(eta1 , -2) - idSF1) / idSF1;// relative syst down uncertainty
    float idSF1D_2 = quad2(idSF1StatD, idSF1SystD); 
    //second muon
    float idSF2StatD = (tnp_weight_muid_pbpb(eta2 , 2) - idSF2) / idSF2;// relative stat down uncertainty
    float idSF2SystD = (tnp_weight_muid_pbpb(eta2 , -2) - idSF2) / idSF2;// relative syst down uncertainty
    float idSF2D_2 = quad2(idSF2StatD, idSF2SystD); 

    idSF_varied = idSF1D_2 + idSF2D_2;//relative uncertainty squared
  }


  //*******************************************************************************************************
  //Muon trigger scale factors:
  //efficiency for data (product of inefficiencies because single muon trigger
  float eff1 = tnp_weight_trig_pbpb(pt1, eta1, 200);
  float eff2 = tnp_weight_trig_pbpb(pt2, eta2, 200);
  float eff_data = 1 - ((1 - eff1)*(1 - eff2));
  float eff_MC = 1 - ( (1 - tnp_weight_trig_pbpb(pt1, eta1, 300) ) * (1 - tnp_weight_trig_pbpb(pt2, eta2, 300)) );

  float trigSF = eff_data/eff_MC;
  float trigSF_varied = trigSF;

  //upward
  if(idx == 1){
    //first uncertainty
    float eff_data_TnPU1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, -1)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_TnPU2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, -1)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) )) - eff_data)/eff_data;
    float TnPU = quad2( eff_data_TnPU1, eff_data_TnPU2);
    //second uncertainty
    float eff_data_StatU1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, 1)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatU2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, 1)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) ))- eff_data)/eff_data; 
    float StatU = quad2( eff_data_StatU1, eff_data_StatU2); 
    //third uncertainty
    float eff_data_RunU1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, 3)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_RunU2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, 3)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) )) - eff_data)/eff_data; 
    float RunU = quad2( eff_data_RunU1, eff_data_RunU2); 

    trigSF_varied = TnPU + StatU + RunU;//relative uncertainty squared
    std::cout << trigSF_varied << std::endl;
  }
  //downward
  if(idx == -1){
    //first uncertainty
    float eff_data_TnPD1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, -2)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_TnPD2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, -2)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) )) - eff_data)/eff_data;
    float TnPD = quad2( eff_data_TnPD1, eff_data_TnPD2); 
    //second uncertainty
    float eff_data_StatD1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, 2)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatD2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, 2)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) ))- eff_data)/eff_data; 
    float StatD = quad2( eff_data_StatD1, eff_data_StatD2); 
    //third uncertainty
    float eff_data_RunD1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, 4)/tnp_weight_trig_pbpb(pt1, eta1, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_RunD2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2, 4)/tnp_weight_trig_pbpb(pt2, eta2, 0) ) )) - eff_data)/eff_data; 
    float RunD = quad2( eff_data_RunD1, eff_data_RunD2); 

    trigSF_varied = TnPD + StatD + RunD;//relative uncertainty squared
  }

  float SF = recoSF * idSF * trigSF;
  std::cout << SF << std::endl;
 
  //check if we are doing a variation
  std::cout << idSF_varied << " " << trigSF_varied << std::endl;
  if(idx==1) return SF * ( 1 + TMath::Sqrt( idSF_varied + trigSF_varied ));//add relative uncertainties on the quantity, then add 1 and multiply by base value
  if(idx==-1) return SF * ( 1 - TMath::Sqrt( idSF_varied + trigSF_varied ));//add relative uncertainties on the quantity, then add 1 and multiply by base value

  return SF;
}

#endif

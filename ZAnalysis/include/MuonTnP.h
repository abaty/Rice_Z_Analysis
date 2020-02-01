#ifndef MUONTNP
#define MUONTNP

#include "tnp_weight.h"
#include "TMath.h"

class MuonTnP{

public:
 float getZSF(float pt1, float eta1, float pt2, float eta2, int cent, int idx, bool doCorrelations = true);
 float getMuonSF(float pt1, float eta1, int cent);
 bool isSameTrkTnPBin(float eta1, float eta2);
 int  findTrkTnPBins(float eta);
 bool isSameIdTnPBin(float eta1, float eta2);
 int  findIdTnPBins(float eta);

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

int MuonTnP::findTrkTnPBins(float eta){
  if( eta >= -2.4 && eta<-2.1) return 0;
  if( eta >= -2.1 && eta<-1.6) return 1;
  if( eta >= -1.6 && eta<-1.2) return 2;
  if( eta >= -1.2 && eta<-0.9) return 3;
  if( eta >= -0.9 && eta<0)    return 4;
  if( eta >= 0 && eta<0.9)     return 5;
  if( eta >= 0.9 && eta<1.2)   return 6;
  if( eta >= 1.2 && eta<1.6)   return 7;
  if( eta >= 1.6 && eta<2.1)   return 8;
  if( eta >= 2.1 && eta<2.4)   return 9;
  return -1;
}

bool MuonTnP::isSameTrkTnPBin(float eta1, float eta2){
  if( findTrkTnPBins( eta1 ) == findTrkTnPBins( eta2) ) return true;
  return false;
}

int MuonTnP::findIdTnPBins(float eta){
  if( eta >= -2.4 && eta<-2.1) return 0;
  if( eta >= -2.1 && eta<-1.6) return 1;
  if( eta >= -1.6 && eta<-1.2) return 2;
  if( eta >= -1.2 && eta<-0.9) return 3;
  if( eta >= -0.9 && eta<-0.6)    return 4;
  if( eta >= -0.6 && eta<-0.3)    return 5;
  if( eta >= -0.3 && eta<0)    return 6;
  if( eta >= 0   && eta<0.3)     return 7;
  if( eta >= 0.3 && eta<0.6)     return 8;
  if( eta >= 0.6 && eta<0.9)     return 9;
  if( eta >= 0.9 && eta<1.2)   return 10;
  if( eta >= 1.2 && eta<1.6)   return 11;
  if( eta >= 1.6 && eta<2.1)   return 12;
  if( eta >= 2.1 && eta<2.4)   return 13;
  return -1;
}

bool MuonTnP::isSameIdTnPBin(float eta1, float eta2){
  if( findIdTnPBins( eta1 ) == findIdTnPBins( eta2) ) return true;
  return false;
}

float MuonTnP::getMuonSF(float pt1, float eta1, int cent){
  float trackSF = tnp_weight_glbPFtrk_pbpb(eta1, cent, 0);
  float idSF = tnp_weight_muid_pbpb(eta1, 0);
  float hltSF = tnp_weight_trig_pbpb(pt1, eta1,cent, 0);
  return trackSF * idSF * hltSF;

}

float MuonTnP::getZSF(float pt1, float eta1, float pt2, float eta2, int cent, int idx, bool doCorrelations){
  
  //************************************************************************************
  //Tracking scale factors
  float trackSF1 =  tnp_weight_glbPFtrk_pbpb(eta1, cent, 0);
  float trackSF2 =  tnp_weight_glbPFtrk_pbpb(eta2, cent, 0);
  
  float trackSF = trackSF1 * trackSF2;
  float trackSF_varied = trackSF;
  
  if(idx == 1){
    //first muon
    float trkSF1StatU = (tnp_weight_glbPFtrk_pbpb(eta1, cent , 1) - trackSF1) / trackSF1;// relative stat up uncertainty
    float trkSF1SystU = (tnp_weight_glbPFtrk_pbpb(eta1, cent , -1) - trackSF1) / trackSF1;// relative syst up uncertainty
    float trkSF1U_2 = quad2(trkSF1StatU, trkSF1SystU); 
    //second muon
    float trkSF2StatU = (tnp_weight_glbPFtrk_pbpb(eta2 , cent,  1) - trackSF2) / trackSF2;// relative stat up uncertainty
    float trkSF2SystU = (tnp_weight_glbPFtrk_pbpb(eta2 , cent, -1) - trackSF2) / trackSF2;// relative syst up uncertainty
    float trkSF2U_2 = quad2(trkSF2StatU, trkSF2SystU);
   
    trackSF_varied = trkSF1U_2 + trkSF2U_2;//relative uncertainty squared

    if( isSameTrkTnPBin( eta1, eta2 ) && doCorrelations) trackSF_varied += 2*trkSF1StatU*trkSF2StatU;//add the remaining part to vary in a correlated fashion if they are in the same eta bin
    if( doCorrelations ) trackSF_varied += 2*trkSF1SystU*trkSF2SystU;
  }
  //downward variation
  if(idx == -1){
    //first muon
    float trkSF1StatD = (tnp_weight_glbPFtrk_pbpb(eta1 , cent,  2) - trackSF1) / trackSF1;// relative stat down uncertainty
    float trkSF1SystD = (tnp_weight_glbPFtrk_pbpb(eta1 , cent, -2) - trackSF1) / trackSF1;// relative syst down uncertainty
    float trkSF1D_2 = quad2(trkSF1StatD, trkSF1SystD); 
    //second muon
    float trkSF2StatD = (tnp_weight_glbPFtrk_pbpb(eta2 , cent, 2) - trackSF2) / trackSF2;// relative stat down uncertainty
    float trkSF2SystD = (tnp_weight_glbPFtrk_pbpb(eta2 , cent,-2) - trackSF2) / trackSF2;// relative syst down uncertainty
    float trkSF2D_2 = quad2(trkSF2StatD, trkSF2SystD); 

    trackSF_varied = trkSF1D_2 + trkSF2D_2;//relative uncertainty squared

    if( isSameTrkTnPBin( eta1, eta2 ) && doCorrelations) trackSF_varied += 2*trkSF1StatD*trkSF2StatD;//add the remaining part to vary in a correlated fashion if they are in the same eta bin
    if( doCorrelations ) trackSF_varied += 2*trkSF1SystD*trkSF2SystD;
  }



  //float trackSF_varied = TMath::Power( (2 * 0.006) , 2);//relative uncertainty (0.006) times 2 because this is correlated, squared

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

    if( isSameIdTnPBin( eta1, eta2 ) && doCorrelations) idSF_varied += 2*idSF1StatU*idSF2StatU;//add the remaining part to vary in a correlated fashion if they are in the same eta bin
    if( doCorrelations ) idSF_varied += 2*idSF1SystU*idSF2SystU;
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

    if( isSameIdTnPBin( eta1, eta2 ) && doCorrelations) idSF_varied += 2*idSF1StatD*idSF2StatD;//add the remaining part to vary in a correlated fashion if they are in the same eta bin
    if( doCorrelations ) idSF_varied += 2*idSF1SystD*idSF2SystD;
  }


  //*******************************************************************************************************
  //Muon trigger scale factors:
  //efficiency for data (product of inefficiencies because single muon trigger
  float eff1 = tnp_weight_trig_pbpb(pt1, eta1, cent, 200);
  float eff2 = tnp_weight_trig_pbpb(pt2, eta2, cent, 200);
  float eff_data = 1 - ((1 - eff1)*(1 - eff2));
  float eff_MC = 1 - ( (1 - tnp_weight_trig_pbpb(pt1, eta1, cent,  300) ) * (1 - tnp_weight_trig_pbpb(pt2, eta2, cent, 300)) );

  float trigSF = eff_data/eff_MC;
  float trigSF_varied = trigSF;

  //upward
  if(idx == 1){
    //first uncertainty
    //assuming uncorrelated
    float eff_data_TnPU1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, cent, -1)/tnp_weight_trig_pbpb(pt1, eta1, cent, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_TnPU2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -1)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0) ) )) - eff_data)/eff_data;
    float TnPU = quad2( eff_data_TnPU1, eff_data_TnPU2);

    //second uncertainty
    float eff_data_StatU1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1,cent, 1)/tnp_weight_trig_pbpb(pt1, eta1,cent, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatU2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, 1)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0) ) ))- eff_data)/eff_data; 
    float StatU = quad2( eff_data_StatU1, eff_data_StatU2); 

    //recalculate if we assume correlations
    if(doCorrelations){
      float eff_data_TnPU1Corr = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, cent, -1)/tnp_weight_trig_pbpb(pt1, eta1, cent, 0) ) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -1)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0)) )) - eff_data)/eff_data;
      TnPU = eff_data_TnPU1Corr * eff_data_TnPU1Corr;

      float eff_data_StatU1Corr = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1,cent, 1)/tnp_weight_trig_pbpb(pt1, eta1,cent, 0) ) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -1)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0)) )) - eff_data)/eff_data;
      StatU = eff_data_StatU1Corr * eff_data_StatU1Corr;
    }
    
    trigSF_varied = TnPU + StatU;//relative uncertainty squared
  }
  //downward
  if(idx == -1){
    //first uncertainty
    float eff_data_TnPD1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1,cent, -2)/tnp_weight_trig_pbpb(pt1, eta1,cent, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_TnPD2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -2)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0) ) )) - eff_data)/eff_data;
    float TnPD = quad2( eff_data_TnPD1, eff_data_TnPD2); 
    //second uncertainty
    float eff_data_StatD1 = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1,cent, 2)/tnp_weight_trig_pbpb(pt1, eta1,cent, 0) ) * (1 - eff2) )) - eff_data)/eff_data; 
    float eff_data_StatD2 = ((1 - ( (1 - eff1) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, 2)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0) ) ))- eff_data)/eff_data; 
    float StatD = quad2( eff_data_StatD1, eff_data_StatD2); 
    
    //recalculate if we assume correlations
    if(doCorrelations){
      float eff_data_TnPD1Corr = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1, cent, -2)/tnp_weight_trig_pbpb(pt1, eta1, cent, 0) ) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -2)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0)) )) - eff_data)/eff_data;
      TnPD = eff_data_TnPD1Corr * eff_data_TnPD1Corr;

      float eff_data_StatD1Corr = ((1 - ( (1 - eff1*tnp_weight_trig_pbpb(pt1, eta1,cent, 2)/tnp_weight_trig_pbpb(pt1, eta1,cent, 0) ) * (1 - eff2*tnp_weight_trig_pbpb(pt2, eta2,cent, -2)/tnp_weight_trig_pbpb(pt2, eta2,cent, 0)) )) - eff_data)/eff_data;
      StatD = eff_data_StatD1Corr * eff_data_StatD1Corr;
    }

    trigSF_varied = TnPD + StatD;//relative uncertainty squared
  }

  float SF =  idSF * trigSF * trackSF;
 
  //check if we are doing a variation
  if(idx==1) return SF * ( 1 + TMath::Sqrt( idSF_varied + trigSF_varied + trackSF_varied + 0.02*0.02 ));//add relative uncertainties on the quantity, then add 1 and multiply by base value
  if(idx==-1) return SF * ( 1 - TMath::Sqrt( idSF_varied + trigSF_varied + trackSF_varied + 0.02*0.02));//add relative uncertainties on the quantity, then add 1 and multiply by base value

  return SF;
}

#endif

#ifndef ELECTRONENERGYSCALE
#define ELECTRONENERGYSCALE

#include "TRandom3.h"
#include "TMath.h"
#include <string>

//eta should be the SuperCluster eta (SCEta)

class ElectronEnergyScale{

public:

  ElectronEnergyScale(std::string mode_);
  ~ElectronEnergyScale();

  float correctPt( float pt, float eta, int hiBin);
  float smearPt(   float pt, float eta, int hiBin);
  float scalePt(   float pt, float eta, int hiBin );

  float varyPt( float pt, float eta, int hiBin, int indx);

private:

  TRandom3 * r;
  std::string mode;

};


float ElectronEnergyScale::varyPt( float pt, float eta, int hiBin, int indx){
  if( TMath::Abs(eta) < 1.442 ){
    if(hiBin<20)       return pt*(1 + indx*0.007);
    else if(hiBin<60)  return pt*(1 + indx*0.008);
    else if(hiBin<200) return pt*(1 + indx*0.009);
  }
  else{
    if(hiBin<20)       return pt*(1 + indx*0.023);
    else if(hiBin<60)  return pt*(1 + indx*0.021);
    else if(hiBin<200) return pt*(1 + indx*0.017);
  }
  return pt;
}


float ElectronEnergyScale::smearPt( float pt, float eta, int hiBin){

 float Z = 91.1876;

 if(mode.compare("MC") == 0 ){
   //EB
   if( TMath::Abs(eta) < 1.442 ){
     if( hiBin >=0 && hiBin <20)    return r->Gaus(1.0, 0.904/Z) * pt;
     if( hiBin >= 20 && hiBin < 60) return r->Gaus(1.0, 1.379/Z) * pt;
     if( hiBin >= 60)               return r->Gaus(1.0, 1.786/Z) * pt;
   } else {//EE
     if( hiBin >=0 && hiBin <20)    return r->Gaus(1.0, 3.051/Z) * pt;
     if( hiBin >= 20 && hiBin < 60) return r->Gaus(1.0, 1.214/Z) * pt;
     if( hiBin >= 60)               return r->Gaus(1.0, 3.451/Z) * pt;
    }
  }
 
  //we shouldn't get here
  return pt; 
}

float ElectronEnergyScale::scalePt( float pt, float eta, int hiBin){
  if(mode.compare("data") == 0){
    //EB
    if( TMath::Abs(eta) < 1.442 ){
      if( hiBin >=0 && hiBin <20)    return 0.99 * pt;
      if( hiBin >= 20 && hiBin < 60) return 1.006 * pt;
      if( hiBin >= 60)               return 1.016 * pt;
    } else {//EE
      if( hiBin >=0 && hiBin <20)    return 0.976 * pt;
      if( hiBin >= 20 && hiBin < 60) return 1.015 * pt;
      if( hiBin >= 60)               return 1.052 * pt;
    }
  } else if(mode.compare("MC") == 0 ){
    //EB
    if( TMath::Abs(eta) < 1.442 ){
      if( hiBin >=0 && hiBin <20)    return 0.974 * pt;
      if( hiBin >= 20 && hiBin < 60) return 0.992 * pt;
      if( hiBin >= 60)               return 1.005 * pt;
    } else {//EE
      if( hiBin >=0 && hiBin <20)    return 0.913 * pt;
      if( hiBin >= 20 && hiBin < 60) return 0.952 * pt;
      if( hiBin >= 60)               return 0.992 * pt;
    }
  }

  //we shouldn't get to here...
  return pt;
}

float ElectronEnergyScale::correctPt( float pt, float eta, int hiBin){
  if(mode.compare("data") == 0){
    return scalePt(pt , eta, hiBin);
  } else if(mode.compare("MC") == 0 ){
    return smearPt( scalePt( pt, eta, hiBin ) , eta, hiBin);  
  }

  return 0.0;
}


ElectronEnergyScale::ElectronEnergyScale( std::string mode_ ){
  mode = mode_;
  r = new TRandom3();
}

ElectronEnergyScale::~ElectronEnergyScale(){
  delete r;
}

#endif

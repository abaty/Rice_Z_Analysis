#ifndef MUONTNP
#define MUONTNP

#include "tnp_weight.h"

class MuonTnP{

public:
 float getZSF(float pt1, float eta1, float pt2, float eta2, int idx);

private:
};

float MuonTnP::getZSF(float pt1, float eta1, float pt2, float eta2, int idx){
  
  //Muon reco scale factors:
  float recoSF1 = 1.0;
  float recoSF2 = 1.0;

  float recoSF = recoSF1 * recoSF2;


  //Muon ID scale factors:
  float idSF1 = tnp_weight_muid_pbpb(eta1 , idx);
  float idSF2 = tnp_weight_muid_pbpb(eta2 , idx);

  float idSF = idSF1 * idSF2;


  //Muon trigger scale factors:
  //efficiency for data (product of inefficiencies because single muon trigger
  float eff_data = 1 - ( (1 - tnp_weight_trig_pbpb(pt1, eta1, 200) ) * (1 - tnp_weight_trig_pbpb(pt2, eta2, 200)) );
  float eff_MC = 1 - ( (1 - tnp_weight_trig_pbpb(pt1, eta1, 300) ) * (1 - tnp_weight_trig_pbpb(pt2, eta2, 300)) );

  float trigSF = eff_data/eff_MC;

  float SF = recoSF * idSF * trigSF;
  return SF;
}

#endif

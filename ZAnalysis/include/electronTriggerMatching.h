#ifndef ELECTRONTRIGMATCH
#define ELECTRONTRIGMATCH
#include "TMath.h"
#include <vector>

class ElecTrigObject{
  public:

    unsigned short L1nEGs;
    std::vector< float > * L1egEta = 0;
    std::vector< float > * L1egPhi = 0;
    std::vector< float > * L1egEt = 0;

    std::vector< double > * HLTEta = 0;
    std::vector< double > * HLTPhi = 0;
    std::vector< double > * HLTPt = 0;
};


class ElectronTriggerMatcher{

public:

  bool isL1Matched(float eleEta, float elePhi, ElecTrigObject trig, float Et = 0);
  bool isHLTMatched(float eleEta, float elePhi, ElecTrigObject trig, float Et = 0);
  bool isMatched(float eleEta, float elePhi, ElecTrigObject trig);
  bool isBothMatched( float eleEta1, float elePhi1, float eleEta2, float elePhi2, ElecTrigObject trig);

private:

  float  dR(float phi1, float eta1, float phi2, float eta2);

};

bool ElectronTriggerMatcher::isL1Matched(float eleEta, float elePhi, ElecTrigObject trig, float Et){
  for(auto i = 0; i < trig.L1nEGs; i++){
    //dR matching
    bool isMatched = ( dR( elePhi, eleEta, trig.L1egPhi->at(i), trig.L1egEta->at(i) ) ) < 0.1;
    //energy matching
    bool passesEt = (trig.L1egEt->at(i) >= Et);
    if( isMatched && passesEt ) return true;
  }
  return false;
}

bool ElectronTriggerMatcher::isHLTMatched(float eleEta, float elePhi, ElecTrigObject trig, float Et){
  for(unsigned int i = 0; i < trig.HLTPhi->size(); i++){
    //energy matching
    bool passesEt = (trig.HLTPt->at(i) >= Et);
    if( !passesEt ) continue;
    //dR matching
    bool isMatched = ( dR( elePhi, eleEta, trig.HLTPhi->at(i), trig.HLTEta->at(i) ) ) < 0.1;
    if( isMatched && passesEt ) return true;
  }
  return false;
}

bool ElectronTriggerMatcher::isMatched(float eleEta, float elePhi, ElecTrigObject trig){
  return isL1Matched( eleEta, elePhi, trig) && isHLTMatched( eleEta, elePhi, trig);
}
  
bool ElectronTriggerMatcher::isBothMatched( float eleEta1, float elePhi1, float eleEta2, float elePhi2, ElecTrigObject trig){
  return isMatched(eleEta1, elePhi1, trig) && isMatched(eleEta2, elePhi2, trig);
}

float ElectronTriggerMatcher::dR(float phi1, float eta1, float phi2, float eta2){
  float dPhi = TMath::ACos(TMath::Cos(TMath::Abs(phi1 - phi2)));
  float dEta = TMath::Abs(eta1 - eta2);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}




#endif

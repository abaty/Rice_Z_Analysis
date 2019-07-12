#ifndef ELECTRONSELECTOR
#define ELECTRONSELECTOR
#include "TMath.h"

class ElectronSelector{
  public:

  enum WorkingPoint{
    veto = 0,
    loose = 1,
    medium = 2,
    tight = 3
  };

  bool isGoodElectron(WorkingPoint wp, int hiBin, float eta, float sigIetaIeta, float dEta, float dPhi, int missHits, float HoverEBc, float EoverPInv, float ip3D);
  

  private:

  float sigIetaIetaCut[4][4] = {{0.0147, 0.0135, 0.0116, 0.0104},{0.048, 0.0466, 0.0418, 0.0358},{0.0113, 0.0107, 0.0101, 0.0099},{0.0376, 0.0339, 0.0316, 0.0288}};
  float dEtaCut[4][4] = {{0.0041, 0.0038, 0.0037, 0.0029},{0.0097, 0.0063, 0.0062, 0.0051},{0.0037, 0.0035, 0.0033, 0.0026},{0.0074, 0.0067, 0.0051, 0.0044}};
  float dPhiCut[4][4] = {{0.0853, 0.0376, 0.0224, 0.0206},{0.2348, 0.1186, 0.0373, 0.0266},{0.1280, 0.0327, 0.0210, 0.0170},{0.2085, 0.0838, 0.0384, 0.0266}};
  float HoverEBcCut[4][4] = {{0.2733, 0.1616, 0.1589, 0.1459},{0.1898, 0.1317, 0.1092, 0.0925},{0.1814, 0.1268, 0.0311, 0.0067},{0.1138, 0.0977, 0.0810, 0.0655}};
  float EoverPInvCut[4][4] = {{0.0367, 0.0177, 0.0173, 0.0105},{0.0300, 0.0201, 0.0133, 0.0065},{0.1065, 0.0774, 0.0701, 0.0077},{0.0237, 0.0193, 0.0192, 0.0123}};
  int missHitsCut[4] = {3, 1, 1, 1};
  float IP3DCut = 0.03;


  int getKinematicIndex(bool isBarrel, bool isCentral);
};

int ElectronSelector::getKinematicIndex(bool isBarrel, bool isCentral){
  if(isBarrel && isCentral) return 0;
  if(!isBarrel && isCentral) return 1;
  if(isBarrel && !isCentral) return 2;
  if(!isBarrel && !isCentral) return 3;
  else return -1;
}

bool ElectronSelector::isGoodElectron(WorkingPoint wp, int hiBin, float eta, float sigIetaIeta, float dEta, float dPhi, int missHits, float HoverEBc, float EoverPInv, float ip3D){
  bool isBarrel = false;
  if(TMath::Abs(eta)<1.442) isBarrel = true;
  else if(TMath::Abs(eta) >= 1.442 && TMath::Abs(eta)<2.4) isBarrel = false;
  else return false;
  
  bool isCentral = false;
  if(hiBin<60) isCentral = true;
  else if(hiBin<200) isCentral = false;
  else return false;

  int indx = getKinematicIndex(isBarrel, isCentral);

  if(sigIetaIeta > sigIetaIetaCut[indx][wp]) return false;
  if(dEta > dEtaCut[indx][wp]) return false;
  if(dPhi > dPhiCut[indx][wp]) return false;
  if(HoverEBc > HoverEBcCut[indx][wp]) return false;
  if(TMath::Abs(EoverPInv) > EoverPInvCut[indx][wp]) return false;
  if(missHits > missHitsCut[wp]) return false;
  if(ip3D > IP3DCut) return false;


  return true;
}

#endif

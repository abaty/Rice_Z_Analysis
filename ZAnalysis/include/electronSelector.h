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

  bool isGoodElectron(WorkingPoint wp, int hiBin, float eta, float sigIetaIeta, float dEta, float dPhi, int missHits, float HoverE, float EoverPInv, float d0, float dz);
  

  private:

  float sigIetaIetaCut[4][4] = {{0.0164, 0.0161, 0.0161, 0.0155},{0.0481, 0.0479, 0.0477, 0.0450},{0.0140, 0.0117, 0.0111, 0.0110},{0.0460, 0.0447, 0.0440, 0.0422}};
  float dEtaCut[4][4] = {{0.0065, 0.0053, 0.0053, 0.0044},{0.0146, 0.0145, 0.0090, 0.0090},{0.0081, 0.0071, 0.0041, 0.0026},{0.0108, 0.0108, 0.0075, 0.0066}};
  float dPhiCut[4][4] = {{0.0610, 0.0288, 0.0147, 0.0094},{0.1025, 0.0516, 0.0334, 0.0174},{0.0302, 0.0221, 0.0157, 0.0107},{0.0576, 0.0301, 0.0191, 0.0126}};
  float HoverECut[4][4] = {{0.3470, 0.1984, 0.1735, 0.1569},{0.1993, 0.1892, 0.1428, 0.1367},{0.2838, 0.1910, 0.1824, 0.1820},{0.2104, 0.1627, 0.1616, 0.1473}};
  float EoverPInvCut[4][4] = {{0.1157, 0.1129, 0.0098, 0.0047},{0.0129, 0.0115, 0.0047, 0.0040},{0.0589, 0.0405, 0.0064, 0.0063},{0.0754, 0.0281, 0.0086, 0.0037}};
  int missHitsCut[4] = {3, 1, 1, 1};
  float d0Cut[4] = {0.01, 0.02, 0.01, 0.02};
  float dzCut[4] = {0.04, 0.04, 0.04, 0.04};


  int getKinematicIndex(bool isBarrel, bool isCentral);
};

int ElectronSelector::getKinematicIndex(bool isBarrel, bool isCentral){
  if(isBarrel && isCentral) return 0;
  if(!isBarrel && isCentral) return 1;
  if(isBarrel && !isCentral) return 2;
  if(!isBarrel && !isCentral) return 3;
  else return -1;
}

bool ElectronSelector::isGoodElectron(WorkingPoint wp, int hiBin, float eta, float sigIetaIeta, float dEta, float dPhi, int missHits, float HoverE, float EoverPInv, float d0, float dz){
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
  if(HoverE > HoverECut[indx][wp]) return false;
  if(TMath::Abs(EoverPInv) > EoverPInvCut[indx][wp]) return false;
  if(missHits > missHitsCut[wp]) return false;
  if(TMath::Abs( d0 ) > d0Cut[wp]) return false;
  if(TMath::Abs( dz ) > dzCut[wp]) return false;


  return true;
}

#endif

#ifndef ZEFFICIENCY
#define ZEFFICIENCY

#include <iostream>
#include <string>
#include "TFile.h"
#include "TEfficiency.h"

class ZEfficiency{

  public:

  ZEfficiency(std::string inputFile, bool isMC_ = false, int variationIndex_ = 0);
  ~ZEfficiency();
  double getEfficiency( double y, double pt, double hiBin);

  private:

  int variationIndex;
  bool isMC;
  TFile * f;
  TEfficiency * e[11];

};

double ZEfficiency::getEfficiency(double y, double pt, double hiBin){
  double originalPt = pt;
  if(pt>=200){
    std::cout << "Very High-Pt Z (>200) detected. Pt = " << originalPt << " I am pretending it's pt is 199.9 for efficiency purposes!" << std::endl;
    pt = 199.9;
  }

  int indx = 0;
  if(hiBin>-1 && hiBin < 10) indx = 0;
  else if(hiBin >= 10 && hiBin < 20) indx = 1;
  else indx = hiBin/20+1;

  int bin = e[indx]->FindFixBin(y, pt);
  float efficiency =  e[indx]->GetEfficiency(bin);

  if(efficiency >=0 && efficiency <=1) return efficiency;
  else{
    std::cout << "efficiency not in the range [0,1], returning 1!" << std::endl;
    std::cout << "Rapidity: " << y << " Pt: " << originalPt << std::endl;
    return 1;
  }
}

ZEfficiency::ZEfficiency(std::string inputFile, bool isMC_, int variationIndex_ ){
  isMC = isMC_;
  variationIndex = variationIndex_;

  std::string suffix = "";
  if(isMC) suffix = "_noSF";
  else if(variationIndex == 1) suffix = "_U";
  else if(variationIndex == -1) suffix = "_D";

  f = TFile::Open(inputFile.c_str(),"read");
  e[0] = (TEfficiency*) f->Get(Form("eff%s_0_5",suffix.c_str()));
  e[1] = (TEfficiency*) f->Get(Form("eff%s_5_10",suffix.c_str()));
  for(int i = 2; i<11; i++){
    e[i] = (TEfficiency*) f->Get(Form("eff%s_%d_%d",suffix.c_str(),10*(i-1),10*i));
  }
}

ZEfficiency::~ZEfficiency(){
  f->Close();
}


#endif

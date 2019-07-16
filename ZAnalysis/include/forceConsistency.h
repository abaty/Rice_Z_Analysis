#ifndef FORCECONSISTENCY
#define FORCECONSISTENCY

#include "TH1.h"
#include "TH2D.h"
#include <iostream>

void forceConsistency(TH2D * pass, TH2D * net){

  for(int i = 0; i< pass->GetSize(); i++){
    if(pass->GetBinContent(i) > net->GetBinContent(i)){
      std::cout << "Warning, fixing an inconsistent bin!" << pass->GetBinContent(i)/net->GetBinContent(i) << std::endl;
      pass->SetBinContent(i,  net->GetBinContent(i));
    }
  }
}

void forceConsistency(TH1D * pass, TH1D * net){

  for(int i = 0; i< pass->GetSize(); i++){
    if(pass->GetBinContent(i) > net->GetBinContent(i)){
      std::cout << "Warning, fixing an inconsistent bin!" << pass->GetBinContent(i)/net->GetBinContent(i) << std::endl;
      pass->SetBinContent(i,  net->GetBinContent(i));
    }
  }
}

#endif

#ifndef HISTHELPER
#define HISTHELPER

#include "TH1D.h"
#include <string>
#include <vector>

class HistNameHelper{
  public:

  HistNameHelper();
  void makeDifferential(TH1D * h);
  ~HistNameHelper();

  std::vector< std::string > name;

};

HistNameHelper::HistNameHelper(){
  name.push_back("massPeak");
  name.push_back("pT");
  name.push_back("y");
  name.push_back("yield");
}

HistNameHelper::~HistNameHelper(){

}

void HistNameHelper::makeDifferential(TH1D* h){
  for(int j = 1; j<h->GetNbinsX()+1; j++){
    h->SetBinContent(j,h->GetBinContent(j)/h->GetBinWidth(j));
    h->SetBinError(j,h->GetBinError(j)/h->GetBinWidth(j));
  }
  return;
}

#endif

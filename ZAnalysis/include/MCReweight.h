#ifndef MCREWEIGHT
#define MCREWEIGHT

#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>

class MCReweight{
  public:

    MCReweight(std::string file, std::string centFile);
    ~MCReweight();
    double reweightFactor(float vz);    
    double reweightFactorCent(int hiBin);

  private:
    TFile * f;
    TH1D * vzRatio;

    TFile * c;
    TH1D * centFlattening;
};

MCReweight::MCReweight(std::string file, std::string centFile){
  f = TFile::Open(file.c_str(),"read");
  vzRatio = (TH1D*)f->Get("vz_Ratio");

  c = TFile::Open(centFile.c_str(),"read");
  centFlattening = (TH1D*)c->Get("DYcent");
}

MCReweight::~MCReweight(){
  f->Close();
  c->Close();
}

double MCReweight::reweightFactor(float vz){
  if(TMath::Abs(vz) > 20){
    std::cout << "Warning, trying to get the vz reweight factor for |vz|>20 cm.  This region is not supported so I am returning 1!!!" << std::endl;
    return 1;
  }

  return vzRatio->GetBinContent( vzRatio->FindBin(vz) );
}

double MCReweight::reweightFactorCent(int hiBin){
  if(hiBin<0 || hiBin>199){
    std::cout << "Warning, bad hiBin number!" << std::endl;
    return 1;
  }

  return 1.0/( centFlattening->GetBinContent( centFlattening->FindBin(hiBin) ) * 200.0);//200 is just to put the weight close to 1

}

#endif

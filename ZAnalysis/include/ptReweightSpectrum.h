#ifndef PTREWEIGHT
#define PTREWEIGHT
#include "TFile.h"
#include "TH1D.h"
#include <string>


class PtReweightSpectrum{
  public:

  PtReweightSpectrum(std::string inputFile);
  ~PtReweightSpectrum();

  double getReweightFactorElectron( float pt , int indx  = 0);
  double getReweightFactorMuon( float pt , int indx  = 0);

  private:

  TFile * f;
  TH1D * eHists[3];
  TH1D * muHists[3];
};

double PtReweightSpectrum::getReweightFactorElectron( float pt, int indx){
  if(indx==-1) indx = 2;

  if(pt>=200) pt = 199;

  return eHists[indx]->GetBinContent(eHists[indx]->FindBin(pt));
}

double PtReweightSpectrum::getReweightFactorMuon( float pt, int indx){
  if(indx==-1) indx = 2;

  if(pt>=200) pt = 199;

  return muHists[indx]->GetBinContent(muHists[indx]->FindBin(pt));
}

PtReweightSpectrum::PtReweightSpectrum(std::string inputFile){
  f = TFile::Open(inputFile.c_str(),"read");
  
  eHists[0] = (TH1D*) f->Get("ptWeightElectrons");
  eHists[1] = (TH1D*) f->Get("ptWeightElectronsU");
  eHists[2] = (TH1D*) f->Get("ptWeightElectronsD");
  muHists[0] = (TH1D*) f->Get("ptWeightMuons");
  muHists[1] = (TH1D*) f->Get("ptWeightMuonsU");
  muHists[2] = (TH1D*) f->Get("ptWeightMuonsD");
}

PtReweightSpectrum::~PtReweightSpectrum(){
  f->Close();
} 

#endif

#ifndef PTCORRECTOR
#define PTCORRECTOR

#include <iostream>
#include <string>
#include "TFile.h"
#include "TEfficiency.h"
#include "TH1D.h"

class ptCorrector{

  public:

  ptCorrector(std::string inputFileMu, std::string inputFileE);
  ~ptCorrector();

  TFile * fmu, *fe;
  TH1D * correction[3];
  TH1D * uncertainty[3];

};


ptCorrector::ptCorrector(std::string inputFileMu, std::string inputFileE ){
  fmu = TFile::Open(inputFileMu.c_str(),"read");
  fe = TFile::Open(inputFileE.c_str(),"read");
  correction[0] = (TH1D*)fmu->Get("recoEff_pt_pass_forReso_Ratio_Reco_0_90");
  correction[1] = (TH1D*)fe->Get("recoEff_pt_pass_forReso_Ratio_Reco_0_90");
  correction[2] = (TH1D*)fmu->Get("recoEff_pt_pass_forReso_Ratio_Reco21_0_90");
  uncertainty[0] = (TH1D*)fmu->Get("recoEff_pt_pass_forReso_Ratio_SmearedToNominal_0_90");
  uncertainty[1] = (TH1D*)fe->Get("recoEff_pt_pass_forReso_Ratio_SmearedToNominal_0_90");
  uncertainty[2] = (TH1D*)fmu->Get("recoEff_pt_pass_forReso_Ratio_SmearedToNominal21_0_90");

  for(int i = 0; i<3; i++){
    for(int j = 0; j<correction[i]->GetSize(); j++){
      correction[i]->SetBinError(j,0);
      uncertainty[i]->SetBinError(j,0);
    }
  }
}

ptCorrector::~ptCorrector(){
  fmu->Close();
  fe->Close();
}


#endif

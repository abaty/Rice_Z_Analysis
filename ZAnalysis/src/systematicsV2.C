#include "include/HistNameHelper.h"
#include "include/MCReweight.h"
#include "include/electronEnergyScale.h"
#include "include/electronSelector.h"
#include "include/electronTriggerMatching.h"
#include "include/centralityBin.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/Timer.h"
#include "include/ZEfficiency.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
TH1D * convertProfileToHistogram(TProfile * input, std::string title){
  int numberOfBins = input->GetXaxis()->GetNbins();
  float lowerBound = input->GetXaxis()->GetBinLowEdge(1);
  float upperBound = input->GetXaxis()->GetBinUpEdge(numberOfBins);
  TH1D * target = new TH1D(title.c_str(),title.c_str(),numberOfBins,lowerBound,upperBound);
  for(int i = 0; i<numberOfBins+2; i++){
    target->SetBinContent(i, input->GetBinContent(i));
    target->SetBinError(i, input->GetBinError(i));
  }
  return target;
}
void sqrtHist(TH1D * h){
  int numberOfBins = h->GetXaxis()->GetNbins();
  for(int i = 0; i<numberOfBins+2; i++){
    h->SetBinContent(i, TMath::Sqrt(h->GetBinContent(i)));
    h->SetBinError(i, 0.5/TMath::Sqrt(h->GetBinContent(i)) ); //   uncertainty on sqrt(x) is: 1/2 * 1/sqrt(x)
  }
}

//void systematicsV2(std::string file, std::string hiBin1, std::string hiBin2, std::string outputTag, bool isEE){
void systematicsV2(std::string file, std::string hiBin1, std::string hiBin2, bool isEE){
  
//  CentralityTool c = CentralityTool();
//  const int nBins = c.getNCentBins();

  HistNameHelper h = HistNameHelper();
 

  TH1D * v2MuMu_Num[3][3] = {0}; 
  TH1D * v2MuMu_Denom[3][3] = {0}; 
  TH1D * v2MuMu_Q1Mid[3][3] = {0}; 
  TH1D * v2MuMu_Q2Mid[3][3] = {0}; 
 
  TProfile * p_v2MuMu_Num[3][3]; 
  TProfile * p_v2MuMu_Denom[3][3]; 
  TProfile * p_v2MuMu_Q1Mid[3][3]; 
  TProfile * p_v2MuMu_Q2Mid[3][3]; 
 
  TFile * ZFile[3];
  ZFile[0] = TFile::Open(file.c_str(),"read");
  ZFile[1] = TFile::Open(hiBin1.c_str(),"read");
  ZFile[2] = TFile::Open(hiBin2.c_str(),"read");
 
  TProfile * pEffCorr[3][3];

  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      if(i>0 && j>0) continue;
      p_v2MuMu_Num[i][j] =   (TProfile*)ZFile[i]->Get(Form("v2NumVsCent%s",h.variationName.at(j).c_str()));
      pEffCorr[i][j] = (TProfile*)ZFile[i]->Get(Form("v2AvgEffVsCent%s",h.variationName.at(j).c_str()));
      v2MuMu_Num[i][j] = convertProfileToHistogram(p_v2MuMu_Num[i][j], "h_v2NumVsCent"+h.variationName.at(j));

      //divide by average efficiency correction to effectively make the efficiency weights all normalized to 1
      for(int k = 0; k<v2MuMu_Num[i][j]->GetSize(); k++){
        v2MuMu_Num[i][j]->SetBinContent(k, v2MuMu_Num[i][j]->GetBinContent(k) / pEffCorr[i][j]->GetBinContent(k) );
        v2MuMu_Num[i][j]->SetBinError(k,   v2MuMu_Num[i][j]->GetBinError(k) / pEffCorr[i][j]->GetBinContent(k) );
      }

      if(j==0){
        p_v2MuMu_Denom[i][j] = (TProfile*)ZFile[i]->Get("v2DenomVsCent");
        v2MuMu_Denom[i][j] = convertProfileToHistogram(p_v2MuMu_Denom[i][j], "h2_v2DenomVsCent"+h.variationName.at(j));
        p_v2MuMu_Q1Mid[i][j] = (TProfile*)ZFile[i]->Get("v2Q1MidVsCent");
        v2MuMu_Q1Mid[i][j] = convertProfileToHistogram(p_v2MuMu_Q1Mid[i][j], "h2_v2Q1MidVsCent"+h.variationName.at(j));
        p_v2MuMu_Q2Mid[i][j] = (TProfile*)ZFile[i]->Get("v2Q2MidVsCent");
        v2MuMu_Q2Mid[i][j] = convertProfileToHistogram(p_v2MuMu_Q2Mid[i][j], "h2_v2Q2MidVsCent"+h.variationName.at(j));
      }
    }
    v2MuMu_Denom[i][0]->Multiply(v2MuMu_Q1Mid[i][0]);
    v2MuMu_Denom[i][0]->Divide(v2MuMu_Q2Mid[i][0]);
    sqrtHist(v2MuMu_Denom[i][0]);
  }
  
  TH1D * v2[3][3];
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      if(i>0 && j>0) continue;
      v2[i][j]= (TH1D*)v2MuMu_Num[i][j]->Clone(Form("v2_MuMu%s_hiBin%d",h.variationName.at(j).c_str(),i));
      v2[i][j]->Divide(v2MuMu_Denom[i][0]);
    }
  }


  /*


 
  TFile * inFile[3];
  inFile[0] = TFile::Open(file.c_str(),"read");
  inFile[1] = TFile::Open(hiBin1.c_str(),"read");
  inFile[2] = TFile::Open(hiBin2.c_str(),"read");
*/
  TH1D * nominal = (TH1D*) v2[0][0]->Clone("v2NumVsCent");
  TH1D * nominalvar1 = (TH1D*) v2[0][1]->Clone("v2NumVsCent_var1");
  TH1D * nominalvar2 = (TH1D*) v2[0][2]->Clone("v2NumVsCent_var2");
  TH1D * nominalhiBin1 = (TH1D*) v2[1][0]->Clone("v2NumVsCent");
  TH1D * nominalhiBin2 = (TH1D*) v2[2][0]->Clone("v2NumVsCent");

  nominalvar1->Add(nominal,-1);
  nominalvar2->Add(nominal,-1);
  nominalhiBin1->Add(nominal,-1);
  nominalhiBin2->Add(nominal,-1);

  h.absAndSymmetrize(nominalvar1, nominalvar2);
  h.absAndSymmetrize(nominalhiBin1, nominalhiBin2);
 
  TH1D * totalError = (TH1D*) nominalvar1->Clone("totalError");
  h.addInQuadrature2(totalError,nominalhiBin1);

  TFile * out = TFile::Open(Form("systematics_v2_isEE%d.root",(int)isEE),"recreate");
  nominalvar1->SetName("effError");
  nominalhiBin1->SetName("hfError");
  nominal->SetName("nominalV2Result");
  nominal->Write();
  nominalvar1->Write();
  nominalhiBin1->Write();
  totalError->Write();

  out->Close();
 
  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: systematics.exe <nominalFile> <hiBin1 File> <hiBin2 File> <isEE>" << std::endl;
    return 1;
  }  

  std::string file  =     argv[1];
  std::string hiBin1    = argv[2];
  std::string hiBin2    = argv[3];
  //std::string outputTag = argv[4];
  bool isEE = (bool)std::atoi(argv[4]);
  systematicsV2(file, hiBin1, hiBin2, isEE);
  //systematicsV2(file, hiBin1, hiBin2, outputTag, isEE);
  return 0; 
}

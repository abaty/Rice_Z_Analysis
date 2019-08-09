#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void dataMCCompare(std::string Zee, std::string ZeeMC, bool isMu){
  //Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS_EE[nBins]; 
  TH1D * massPeakOS_EEMC[nBins]; 
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZeeMCFile = TFile::Open(ZeeMC.c_str(),"read");
    
  for(int i = 0; i<nBins; i++){
    massPeakOS_EE[i] = (TH1D*)ZeeFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_EEMC[i] = (TH1D*)ZeeMCFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    float norm =  massPeakOS_EE[i]->Integral();
    massPeakOS_EE[i]->Scale(1.0/norm);
    norm =  massPeakOS_EEMC[i]->Integral();
    massPeakOS_EEMC[i]->Scale(1.0/norm);


    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetBottomMargin(0.2);
    massPeakOS_EE[i]->GetYaxis()->SetRangeUser(0, massPeakOS_EEMC[i]->GetMaximum()*1.6);  
    massPeakOS_EE[i]->GetXaxis()->SetTitle("m_{ll} (GeV)");
    massPeakOS_EE[i]->GetYaxis()->SetTitle("Counts");
    massPeakOS_EE[i]->GetXaxis()->CenterTitle();
    massPeakOS_EE[i]->GetYaxis()->CenterTitle();

    massPeakOS_EE[i]->SetStats(0);
    
    massPeakOS_EE[i]->SetLineColor(kBlack);
    massPeakOS_EE[i]->SetFillColor(kOrange-2);
    massPeakOS_EE[i]->Draw("HIST");
    massPeakOS_EEMC[i]->SetLineColor(kBlack);
    massPeakOS_EEMC[i]->SetFillColor(kBlack);
    massPeakOS_EEMC[i]->SetFillStyle(3345);
    massPeakOS_EEMC[i]->Draw("HIST same");
    
    TLegend * leg = new TLegend(0.28,0.72,0.63,0.87);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if(!isMu){
      leg->AddEntry(massPeakOS_EE[i],"Z #rightarrow e^{+}e^{-} Data","f");
      leg->AddEntry(massPeakOS_EEMC[i],"Z #rightarrow e^{+}e^{-} MC","f");
    }
    else{
      leg->AddEntry(massPeakOS_EE[i],"Z #rightarrow #mu^{+}#mu^{-} Data","f");
      leg->AddEntry(massPeakOS_EEMC[i],"Z #rightarrow #mu^{+}#mu^{-} MC","f");
    }

    leg->Draw("same");
    
    TLegend * leg2 = new TLegend(0.65,0.72+0.15*0.25,0.8,0.87);
    leg2->SetBorderSize(0);
    leg2->AddEntry((TObject*)0,"2018 PbPb","");
    leg2->AddEntry((TObject*)0,Form("%d-%d %%",c.getCentBinLow(i), c.getCentBinHigh(i)),"");
    leg2->AddEntry((TObject*)0,"p_{T}^{l} > 20 GeV","");
    leg2->Draw("same");


    c1->SaveAs(Form("plots/dataVsMC/massPeak_isMu%d_%d_%d.png", (int)isMu ,c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/dataVsMC/massPeak_isMu%d_%d_%d.pdf", (int)isMu ,c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/dataVsMC/massPeak_isMu%d_%d_%d.C", (int)isMu ,c.getCentBinLow(i),c.getCentBinHigh(i)));
  
    delete leg2; 
    delete leg;
    delete c1;
  }
  return;

  
}

int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: dataMCCompare.exe <data file> <MC file> <isMu>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string ZeeMC = argv[2];
  bool isMu = (bool)std::atoi(argv[3]);
   
  dataMCCompare(Zee, ZeeMC, isMu);
  return 0; 
}

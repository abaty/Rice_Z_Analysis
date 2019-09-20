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

void systematics(std::string file, std::string hiBin1, std::string hiBin2, std::string effTable, std::string outputTag, bool isMu21, bool isEE){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  HistNameHelper h = HistNameHelper();
  //Settings s = Settings();
  
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();
  
  TFile * inFile[3];
  inFile[0] = TFile::Open(file.c_str(),"read");
  inFile[1] = TFile::Open(hiBin1.c_str(),"read");
  inFile[2] = TFile::Open(hiBin2.c_str(),"read");

  TH1D * result[nBins][4][5][3];
  TH1D * backgroundYields[nBins][4][3];//nBin cent bins, 4 observables, 3 backgrounds

  //get a lot of stuff from the variations
  for(int f = 0; f< 3; f++){
    for(int i = 0; i<nBins; i++){
      for(int j = 1; j<4; j++){//skip j==0 (masses)
        for(int k = 0; k<5; k++){
          if(f>0 && k>0) continue;//skip variations of centrality variations
  
          result[i][j][k][f] = (TH1D*)inFile[f]->Get(Form("%sOS_minusAll%s_%d_%d",h.name.at(j).c_str(),h.variationName.at(k).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)));

          if(f==0){
            backgroundYields[i][j][0] = (TH1D*) inFile[f]->Get(Form("%sBkg_Wjet_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
            backgroundYields[i][j][1] = (TH1D*) inFile[f]->Get(Form("%sBkg_ttbar_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
            backgroundYields[i][j][2] = (TH1D*) inFile[f]->Get(Form("%sBkg_tau_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
          }
        }
      }    
    }
  }


  std::cout << "Dividing things now!" << std::endl;
  //lets divide things by the nominal result now, first the centrality stuff
  for(int i = 0; i<nBins; i++){
    for(int j = 1; j<4; j++){
      result[i][j][0][1]->Divide(result[i][j][0][0]);//cent up, include Nmb scaling here
      result[i][j][0][1]->Scale(1.01);

      result[i][j][0][2]->Divide(result[i][j][0][0]);//cent down, include Nmb scaling here
      result[i][j][0][2]->Scale(0.995);
    }
  }

  std::cout << "Other Variations now!" << std::endl;
  //lets do the other 4 variations that are in the same file
  for(int i = 0; i<nBins; i++){
    for(int j = 1; j<4; j++){
      for(int k = 1; k<5; k++){
        result[i][j][k][0]->Divide(result[i][j][0][0]);
      }     
    }
  }
  
  std::cout << "Backgrounds now!" << std::endl;
  //make the background yields into relative differences
  for(int i = 0; i<nBins; i++){
    for(int j = 1; j<4; j++){
      for(int k = 0; k<3; k++){//3 backgrounds here
        backgroundYields[i][j][k]->Scale(0.2);//we assume they vary by 10\%
        backgroundYields[i][j][k]->Divide(result[i][j][0][0]);
      }     
    }
  }

  std::cout << "MC stat uncertainties now!" << std::endl;
  TH1D * mcStatRelErr[nBins];
  TH1D * mcStatRelErr_pt;
  TH1D * mcStatRelErr_y;
  for(int i = 0; i<nBins; i++){
    mcStatRelErr[i] = (TH1D*) inFile[0]->Get(Form("yieldOS_withEff_RelStatErr_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
  } 
  mcStatRelErr_pt = (TH1D*)inFile[0]->Get("pTOS_withEff_RelStatErr_0_90");
  mcStatRelErr_y = (TH1D*)inFile[0]->Get("yOS_withEff_RelStatErr_0_90");


  //get pt smearing effect
  std::cout << "Pt Smearing now!" << std::endl;
  std::vector< std::string > label21;
  label21.push_back("21");
  label21.push_back("");
  //TH1D * reso[nBins];
  TFile * eff = TFile::Open(effTable.c_str(),"read"); 
  //for(int i = 0; i<nBins; i++){
  //  reso[i] = (TH1D*) eff->Get(Form("recoEff_pt_pass_forReso_Ratio_SmearedToNominal%s_%d_%d", isMu21 ? label21[0].c_str() : label21[1].c_str() ,c.getCentBinLow(i),c.getCentBinHigh(i)));
  //}
  
  TH1D * accept[4];
  accept[1] = (TH1D*)eff->Get(Form("accept%d_totalUncert_pt", (isMu21 || isEE) ? 21 : 24));
  accept[2] = (TH1D*)eff->Get(Form("accept%d_totalUncert_y", (isMu21 || isEE) ? 21 : 24 ));
  accept[3] = (TH1D*)eff->Get(Form("accept%d_totalUncert_yields", (isMu21 || isEE) ? 21 : 24 ));

  std::cout << "Opening Output File." << std::endl;
  TFile * output = TFile::Open(Form("systematics_%s_isMu21%d.root",outputTag.c_str(),(int)isMu21),"recreate");

  TH1D * effError[nBins][4];
  TH1D * acoError[nBins][4];
  TH1D * hfUError[nBins][4];
  TH1D * hfDError[nBins][4];
  TH1D * hfError[nBins][4];
  TH1D * ptSmearError[nBins][4];
  TH1D * wError[nBins][4];
  TH1D * ttbarError[nBins][4];
  TH1D * tauError[nBins][4];
  TH1D * mcStatError[nBins][4];
  TH1D * chargeSwapError[nBins][4];
  TH1D * acceptError[nBins][4];
  TH1D * totalError[nBins][4];


  for(int i = 0; i<nBins; i++){
    for(int j = 1; j<4; j++){
      effError[i][j] = (TH1D*)result[i][j][1][0]->Clone(Form("%s_efficiencyError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      h.subtractOneAndAbsAndSymmetrize(effError[i][j], result[i][j][2][0]);
      
      acoError[i][j] = (TH1D*)result[i][j][3][0]->Clone(Form("%s_emError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      h.subtractOneAndAbsAndSymmetrize(acoError[i][j], result[i][j][4][0]);
      
      hfUError[i][j] = (TH1D*)result[i][j][0][1]->Clone(Form("%s_hfUError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      hfDError[i][j] = (TH1D*)result[i][j][0][2]->Clone(Form("%s_hfDError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      hfError[i][j] = (TH1D*)result[i][j][0][1]->Clone(Form("%s_hfError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      h.subtractOneAndAbsAndSymmetrize(hfError[i][j], result[i][j][0][2]);

      if(j==1){
        //ptSmearError[i][j] = (TH1D*) reso[i]->Clone(Form("%s_ptSmearError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
        //h.subtractOneAndAbs(ptSmearError[i][j]);
        ptSmearError[i][j] = (TH1D*) inFile[0]->Get("unfoldingUncert");
        ptSmearError[i][j]->SetName(Form("%s_ptSmearError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      }

      wError[i][j] = (TH1D*) backgroundYields[i][j][0]->Clone(Form("%s_wError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      ttbarError[i][j] = (TH1D*) backgroundYields[i][j][1]->Clone(Form("%s_ttbarError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      tauError[i][j] = (TH1D*) backgroundYields[i][j][2]->Clone(Form("%s_tauError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      
      //charge Swapping
      if(isEE){
        chargeSwapError[i][j] = (TH1D*)backgroundYields[i][j][0]->Clone(Form("%s_chargeSwapError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
        chargeSwapError[i][j]->Reset();
        for(int k = 0; k<chargeSwapError[i][j]->GetSize(); k++){
          chargeSwapError[i][j]->SetBinContent(k,0.005);
        }
      }   

      //acceptance
      acceptError[i][j] = (TH1D*) effError[i][j]->Clone(Form("%s_acceptError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i), c.getCentBinHigh(i)));
      acceptError[i][j]->Reset();
      for(int k = 0; k<acceptError[i][j]->GetSize(); k++){
        acceptError[i][j]->SetBinContent(k, accept[j]->GetBinContent( accept[j]->FindBin( acceptError[i][j]->GetBinCenter(k) ) ) );
      } 

      //MC Stats
      if(j==1 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) mcStatError[i][j] = (TH1D*) mcStatRelErr_pt->Clone(Form("%s_mcStatError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i))); 
      if(j==2 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) mcStatError[i][j] = (TH1D*) mcStatRelErr_y->Clone(Form("%s_mcStatError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i))); 
      if(j==3) mcStatError[i][j] = (TH1D*) mcStatRelErr[i]->Clone(Form("%s_mcStatError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));

      totalError[i][j] = (TH1D*)effError[i][j]->Clone(Form("%s_totalError_%d_%d",h.name.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)));
      
      if(j!=1) h.addInQuadrature3( totalError[i][j], acoError[i][j], hfError[i][j]);
      if(j==1) h.addInQuadrature4( totalError[i][j], acoError[i][j], hfError[i][j], ptSmearError[i][j]);
 
      if(j==1 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) h.addInQuadrature2( totalError[i][j], mcStatError[i][j]); 
      if(j==2 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) h.addInQuadrature2( totalError[i][j], mcStatError[i][j]); 
      if(j==3) h.addInQuadrature2( totalError[i][j], mcStatError[i][j]);
      
      h.addInQuadrature2(totalError[i][j], acceptError[i][j]);

      if(isEE) h.addInQuadrature2( totalError[i][j], chargeSwapError[i][j]);
    }        
  }

  std::cout << "Writing" << std::endl;

  const char * labels[9] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%", "50-70%", "70-90%", "0-90%"};
  int binMap[9] = {0,1,2,3,4,5,15,16,25};

  TH1D * effErrorByCent = new TH1D("effErrorByCent","",9,0,9);
  TH1D * acoErrorByCent = new TH1D("acoErrorByCent","",9,0,9);
  TH1D * hfErrorByCent = new TH1D("hfErrorByCent","",9,0,9);
  TH1D * totalErrorByCent = new TH1D("totalErrorByCent","",9,0,9);
  TH1D * mcStatErrorByCent = new TH1D("mcStatErrorByCent","",9,0,9);
  TH1D * chargeSwapErrorByCent = new TH1D("chargeSwapErrorByCent","",9,0,9); 
  TH1D * acceptErrorByCent = new TH1D("acceptErrorByCent","",9,0,9);
 
  for(int i = 0; i<nBins; i++){
    for(int j = 1; j<4; j++){
      effError[i][j]->Write();
      acoError[i][j]->Write();
      hfError[i][j]->Write();
      hfUError[i][j]->Write();
      hfDError[i][j]->Write();
      if(j==1) ptSmearError[i][j]->Write();
      if(j==1 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) mcStatError[i][j]->Write(); 
      if(j==2 && c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) mcStatError[i][j]->Write(); 
      if(j==3) mcStatError[i][j]->Write();
      wError[i][j]->Write();
      ttbarError[i][j]->Write();
      tauError[i][j]->Write();
      acceptError[i][j]->Write();
      if(isEE) chargeSwapError[i][j]->Write();
      totalError[i][j]->Write();

      int binIdx = -1;
      if(j==3){
        for(int k = 0; k<9; k++){
          if( binMap[k] == i ) binIdx = k;
        }
        if( binIdx < 0 ) continue;
        
        effErrorByCent->SetBinContent( effErrorByCent->FindBin(binIdx) , effError[i][j]->GetBinContent(1));
        acoErrorByCent->SetBinContent( acoErrorByCent->FindBin(binIdx) , acoError[i][j]->GetBinContent(1));
        hfErrorByCent->SetBinContent( hfErrorByCent->FindBin(binIdx) , hfError[i][j]->GetBinContent(1));
        mcStatErrorByCent->SetBinContent( mcStatErrorByCent->FindBin(binIdx), mcStatError[i][j]->GetBinContent(1));
        acceptErrorByCent->SetBinContent( acceptErrorByCent->FindBin(binIdx), acceptError[i][j]->GetBinContent(1));
        if(isEE) chargeSwapErrorByCent->SetBinContent(chargeSwapErrorByCent->FindBin(binIdx), chargeSwapError[i][j]->GetBinContent(1));
        totalErrorByCent->SetBinContent( totalErrorByCent->FindBin(binIdx) , totalError[i][j]->GetBinContent(1));
      }
 
      TCanvas * c1 = new TCanvas("","",800,800);
      if(j==1) totalError[i][j]->GetXaxis()->SetTitle("p_{T}");
      if(j==2) totalError[i][j]->GetXaxis()->SetTitle("y");
      
      totalError[i][j]->GetYaxis()->SetTitle("Relative Sys. Uncertainty");
      totalError[i][j]->GetYaxis()->SetRangeUser(0,1.7*totalError[i][j]->GetMaximum());
      totalError[i][j]->SetLineColor(kBlack);
      totalError[i][j]->SetStats(0);
      totalError[i][j]->Draw();
      effError[i][j]->SetLineColor(kRed);
      effError[i][j]->Draw("same");
      acoError[i][j]->SetLineColor(kBlue);
      acoError[i][j]->Draw("same");
      hfError[i][j]->SetLineColor(kGreen);
      hfError[i][j]->Draw("same");
      if(j==1){
        ptSmearError[i][j]->SetLineColor(kViolet);
        ptSmearError[i][j]->Draw("same");
      }

      if( ( (j==1 || j==2) && (c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) ) || j==3 ){
        mcStatError[i][j]->SetLineColor(kOrange);
        mcStatError[i][j]->Draw("same");
      }     

      if(isEE){
        chargeSwapError[i][j]->SetLineColor(kGray);
        chargeSwapError[i][j]->Draw("same");
      }

      acceptError[i][j]->SetLineColor(kRed);
      acceptError[i][j]->SetLineStyle(10);
      acceptError[i][j]->Draw("same");

      TLegend * leg = new TLegend(0.3,0.6,0.8,0.9);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      if(isEE) leg->AddEntry((TObject*)0, "ee |#eta_{l} |<2.1","");
      else if(isMu21) leg->AddEntry((TObject*)0, "#mu#mu |#eta_{l} |<2.1","");
      else leg->AddEntry((TObject*)0, "#mu#mu |#eta_{l} |<2.4","");

      leg->AddEntry(totalError[i][j], "Total Uncertainty","l");
      leg->AddEntry(effError[i][j], "Lepton Reco","l");
      leg->AddEntry(acoError[i][j],"EM Background Cut","l");
      leg->AddEntry(hfError[i][j], "HF Uncert. (N_{MB} + centrality)","l");
      if(j==1) leg->AddEntry(ptSmearError[i][j], "p_{T} unfolding","l");
      if( ( (j==1 || j==2) && (c.getCentBinLow(i)==0 && c.getCentBinHigh(i)==90) ) || j==3 ){
        leg->AddEntry(mcStatError[i][j],"MC stat uncertainty","l");
      }
      if(isEE) leg->AddEntry(chargeSwapError[i][j],"Charge Swapping","l");
      leg->AddEntry(acceptError[i][j], "Acceptance","l");
      leg->AddEntry((TObject*)0,"MC Background uncert negligible","");
      if(j==3) leg->AddEntry((TObject*)0,"Glauber Uncert not shown","");
      leg->Draw("same");
   
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d.png",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d.pdf",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d.C",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
      if(j==1) c1->SetLogx();
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d_log.png",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d_log.pdf",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
      c1->SaveAs(Form("plots/systematics/%s_syst_isMu21%d_isMu%d_%d_%d_log.C",h.name.at(j).c_str(),(int)isMu21,(int)!isEE, c.getCentBinLow(i),c.getCentBinHigh(i)));
     

      delete c1;
      delete leg;

    }
  }
  effErrorByCent->Print("All");
  effErrorByCent->Write();   
  acoErrorByCent->Write();   
  hfErrorByCent->Write();  
  acceptErrorByCent->Write(); 
  if(isEE) chargeSwapErrorByCent->Write(); 
  totalErrorByCent->Write();

  TCanvas * c1 = new TCanvas("","",800,800);
  totalErrorByCent->GetXaxis()->SetTitle("Centrality");
  totalErrorByCent->GetYaxis()->SetTitle("Relative Sys. Uncertainty");
  totalErrorByCent->GetYaxis()->SetRangeUser(0,0.15);
  totalErrorByCent->SetLineColor(kBlack);
  totalErrorByCent->SetStats(0);
  for(int i = 1; i<10; i++){
    totalErrorByCent->GetXaxis()->SetBinLabel(i, labels[i-1]);
    totalErrorByCent->GetXaxis()->ChangeLabel(i,45);
  }
  totalErrorByCent->Draw();
  effErrorByCent->SetLineColor(kRed);
  effErrorByCent->Draw("same");
  acoErrorByCent->SetLineColor(kBlue);
  acoErrorByCent->Draw("same");
  hfErrorByCent->SetLineColor(kGreen);
  hfErrorByCent->Draw("same");
  mcStatErrorByCent->SetLineColor(kOrange);
  mcStatErrorByCent->Draw("same");
  chargeSwapErrorByCent->SetLineColor(kGray);
  if(isEE) chargeSwapErrorByCent->Draw("same"); 
  acceptErrorByCent->SetLineColor(kRed);
  acceptErrorByCent->SetLineStyle(10);
  acceptErrorByCent->Draw("same");
 
  TLegend * leg = new TLegend(0.3,0.6,0.8,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(isEE) leg->AddEntry((TObject*)0, "ee |#eta_{l} |<2.1","");
  else if(isMu21) leg->AddEntry((TObject*)0, "#mu#mu |#eta_{l} |<2.1","");
  else leg->AddEntry((TObject*)0, "#mu#mu |#eta_{l} |<2.4","");

  leg->AddEntry(totalErrorByCent, "Total Uncertainty","l");
  leg->AddEntry(effErrorByCent, "Lepton Reco","l");
  leg->AddEntry(acoErrorByCent,"EM Background Cut","l");
  leg->AddEntry(hfErrorByCent, "HF Uncert. (N_{MB} + centrality)","l");
  leg->AddEntry(mcStatErrorByCent,"MC stat uncertainty","l");
  if(isEE) leg->AddEntry(chargeSwapErrorByCent,"Charge Swapping","l");
  leg->AddEntry(acceptErrorByCent,"Acceptance","l");
  leg->AddEntry((TObject*)0,"MC Background uncert negligible","");
  leg->AddEntry((TObject*)0,"Glauber Uncert not shown","");
  leg->Draw("same");
  c1->SaveAs(Form("plots/systematics/%sByCent_syst_isMu21%d_isMu%d.png",h.name.at(3).c_str(),(int)isMu21,(int)!isEE));
  c1->SaveAs(Form("plots/systematics/%sByCent_syst_isMu21%d_isMu%d.pdf",h.name.at(3).c_str(),(int)isMu21,(int)!isEE));
  c1->SaveAs(Form("plots/systematics/%sByCent_syst_isMu21%d_isMu%d.C",h.name.at(3).c_str(),(int)isMu21,(int)!isEE));
  
 
  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 8)
  {
    std::cout << "Usage: systematics.exe <nominalFile> <hiBin1 File> <hiBin2 File> <effTable> <tag> <isMu21> <isEE>" << std::endl;
    return 1;
  }  

  std::string file  =     argv[1];
  std::string hiBin1    = argv[2];
  std::string hiBin2    = argv[3];
  std::string effTable  = argv[4];
  std::string outputTag = argv[5];
  bool isMu21 = (bool)std::atoi(argv[6]);
  bool isEE = (bool)std::atoi(argv[7]);
  systematics(file, hiBin1, hiBin2, effTable, outputTag, isMu21, isEE);
  return 0; 
}

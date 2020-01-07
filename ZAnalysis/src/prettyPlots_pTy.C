//#include "include/ptCorrector.h"
#include "TExec.h"
#include "include/CMS_lumi.C"
#include "include/combinePoints.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/HistNameHelper.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void prettyPlots(std::string Zee, std::string Zmumu21, std::string Zmumu24, std::string ZeeSyst, std::string Zmumu21Syst, std::string Zmumu24Syst, bool doAccept){
  Settings s = Settings();
  CombinePoints cp = CombinePoints();
  HistNameHelper helper = HistNameHelper();

  //ptCorrector ptCorr = ptCorrector("resources/Z2mumu_Efficiencies.root","resources/Z2ee_EfficiencyMC_0.root");

  TH1D * acceptE[2];
  TH1D * acceptMu21[2];
  TH1D * acceptMu24[2];
  if(doAccept){
    TFile * f_acceptE = TFile::Open("resources/Z2ee_EfficiencyMC_0.root");
    acceptE[0] = (TH1D*) f_acceptE->Get("accept21_pt_ratio_0");
    acceptE[1] = (TH1D*) f_acceptE->Get("accept21_y_ratio_0");
    TFile * f_acceptMu = TFile::Open("resources/Z2mumu_Efficiencies.root");
    acceptMu21[0] = (TH1D*) f_acceptMu->Get("accept21_pt_ratio_0");
    acceptMu21[1] = (TH1D*) f_acceptMu->Get("accept21_y_ratio_0");
    acceptMu24[0] = (TH1D*) f_acceptMu->Get("accept24_pt_ratio_0");
    acceptMu24[1] = (TH1D*) f_acceptMu->Get("accept24_y_ratio_0");

    for(int i = 0; i<2; i++){
      for(int j = 0; j<acceptE[i]->GetSize(); j++)  acceptE[i]->SetBinError(j,0);
      for(int j = 0; j<acceptMu21[i]->GetSize(); j++)  acceptMu21[i]->SetBinError(j,0);
      for(int j = 0; j<acceptMu24[i]->GetSize(); j++)  acceptMu24[i]->SetBinError(j,0);
    }
  }

  //CentralityTool c = CentralityTool();
  //const int nBins = c.getNCentBins();

  TFile * theory = TFile::Open("resources/Z2ee_EfficiencyMC_0.root","read");
  TH1D * EPPS16Pt = (TH1D*)theory->Get("theoryEPPS16_pt_Band");
  TH1D * EPPS16Rap = (TH1D*)theory->Get("theoryEPPS16_rap_Band");
  TH1D * nCTEQ15Pt = (TH1D*)theory->Get("theorynCTEQ15_pt_Band");
  TH1D * nCTEQ15Rap = (TH1D*)theory->Get("theorynCTEQ15_rap_Band");
  EPPS16Pt->SetFillColor(kGray+1);
  EPPS16Pt->SetLineColor(kGray+1);
  EPPS16Pt->SetMarkerColor(kGray+1);
  EPPS16Rap->SetFillColor(kGray+1);
  EPPS16Rap->SetLineColor(kGray+1);
  EPPS16Rap->SetMarkerColor(kGray+1);
  nCTEQ15Pt->SetFillColor(kGreen+2);
  nCTEQ15Pt->SetFillStyle(3244);
  nCTEQ15Pt->SetLineColor(kGreen+2);
  nCTEQ15Pt->SetMarkerSize(0);
  nCTEQ15Rap->SetFillColor(kGreen+2);
  nCTEQ15Rap->SetFillStyle(3244);
  nCTEQ15Rap->SetLineColor(kGreen+2);
  nCTEQ15Rap->SetMarkerSize(0);
  
  //gStyle->SetErrorX(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //mumu channel first
  TFile * mu24 = TFile::Open(Zmumu24.c_str(),"read");
  TH1D * y_mu24 = (TH1D*)mu24->Get("yOS_minusAll_0_90");
  TH1D * pt_mu24 = (TH1D*)mu24->Get("unfoldedDist");

  //pt_mu24->Divide( ptCorr.correction[0] );
 
  TFile * mu21 = TFile::Open(Zmumu21.c_str(),"read");
  TH1D * y_mu21 = (TH1D*)mu21->Get("yOS_minusAll_0_90");
  TH1D * pt_mu21 = (TH1D*)mu21->Get("unfoldedDist");
  
  //pt_mu21->Divide( ptCorr.correction[2] );
  
  TFile * e = TFile::Open(Zee.c_str(),"read");
  TH1D * y_e = (TH1D*)e->Get("yOS_minusAll_0_90");
  TH1D * pt_e = (TH1D*)e->Get("unfoldedDist");
  
  //pt_e->Divide( ptCorr.correction[1] );

  //**********************************************Systematics********************************************************
  TFile * systFile[3];
  systFile[0] = TFile::Open(ZeeSyst.c_str(),"read"); 
  systFile[1] = TFile::Open(Zmumu21Syst.c_str(),"read"); 
  systFile[2] = TFile::Open(Zmumu24Syst.c_str(),"read"); 
  
  TH1D * efficiencyError_0_90[3][3];
  TH1D * emError_0_90[3][3];	
  TH1D * hfError_0_90[3][3];	
  TH1D * ptSmearError_0_90[3][3];
  TH1D * mcStatError_0_90[3][3];	
  TH1D * totalError_0_90[3][3];	
  TH1D * acceptError_0_90[3][3];

  TH1D * combo[3];
  combo[1] = (TH1D*) pt_e->Clone("combined_pT_0_90"); 
  combo[2] = (TH1D*) y_e->Clone("combined_y_0_90"); 
  TH1D * comboSyst[3];
  comboSyst[1] = (TH1D*) pt_e->Clone("combinedSyst_pT_0_90"); 
  comboSyst[2] = (TH1D*) y_e->Clone("combinedSyst_y_0_90"); 

  for(int i = 0; i<3; i++){
    for(int j = 1; j<3; j++){
      std::cout << i << " " << j << std::endl;
      efficiencyError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_efficiencyError_0_90",helper.name.at(j).c_str()));
      emError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_emError_0_90",helper.name.at(j).c_str()));
      hfError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_hfError_0_90",helper.name.at(j).c_str()));
      if(j==1) ptSmearError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_ptSmearError_0_90",helper.name.at(j).c_str()));
      mcStatError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_mcStatError_0_90",helper.name.at(j).c_str()));
      acceptError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_acceptError_0_90",helper.name.at(j).c_str()));
      totalError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_totalError_0_90",helper.name.at(j).c_str()));
    }
  }
  //plots here

  
  //scaling
  y_mu24->Scale(s.netLumi/(s.muLumi * s.Nmb));
  y_mu21->Scale(s.netLumi/(s.muLumi * s.Nmb));
  y_e->Scale(s.netLumi/(s.eLumi * s.Nmb));
  pt_mu24->Scale(s.netLumi/(s.muLumi * s.Nmb));
  pt_mu21->Scale(s.netLumi/(s.muLumi * s.Nmb));
  pt_e->Scale(s.netLumi/(s.eLumi * s.Nmb));

  if(doAccept){
    pt_e->Divide(acceptE[0]);
    pt_mu21->Divide(acceptMu21[0]);
    pt_mu24->Divide(acceptMu24[0]);
    y_e->Divide(acceptE[1]);
    y_mu21->Divide(acceptMu21[1]);
    y_mu24->Divide(acceptMu24[1]);
  }

  //combination
  for(int j = 1; j<3; j++){
    TH1D * tempHistE, *tempHistMu;
    if(j==1){
      tempHistE = pt_e;
      tempHistMu = pt_mu21;
    }
    if(j==2){
      tempHistE = y_e;
      tempHistMu = y_mu21;
    }

    tempHistE->Print("All");
    tempHistMu->Print("All");

    for(int i = 0; i<combo[j]->GetSize()+2; i++){
      float e = tempHistE->GetBinContent(i);
      float eStatErr = tempHistE->GetBinError(i);
      float eEffErr = efficiencyError_0_90[0][j]->GetBinContent(i) * e;
      float eEmErr = emError_0_90[0][j]->GetBinContent(i) * e;
      float eHfErr = hfError_0_90[0][j]->GetBinContent(i) * e;
      float eMCStatErr = mcStatError_0_90[0][j]->GetBinContent(i) * e;
      float ePtSmearErr = 0;
      if(j==1) ePtSmearErr = ptSmearError_0_90[0][j]->GetBinContent(i) * e;
      float eChargeSwapErr = 0.005 * e;
      float eAcceptErr = acceptError_0_90[0][j]->GetBinContent(i) * e;

      
      float mu = tempHistMu->GetBinContent(i);
      float muStatErr = tempHistMu->GetBinError(i);
      float muEffErr = efficiencyError_0_90[1][j]->GetBinContent(i) * mu;
      float muEmErr = emError_0_90[1][j]->GetBinContent(i) * mu;
      float muHfErr = hfError_0_90[1][j]->GetBinContent(i) * mu;
      float muMCStatErr = mcStatError_0_90[1][j]->GetBinContent(i) * mu;
      float muPtSmearErr = 0;
      if(j==1) muPtSmearErr = ptSmearError_0_90[1][j]->GetBinContent(i) * mu;
      float muAcceptErr = acceptError_0_90[1][j]->GetBinContent(i) * mu;

      double scaleFactor = 10000000;     
 
      std::vector< TMatrixD > covariance;
      covariance.push_back( cp.getFullUncorrMatrix(muStatErr*scaleFactor, eStatErr*scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(0.0, eChargeSwapErr*scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(muEffErr*scaleFactor, eEffErr*scaleFactor) );
      covariance.push_back( cp.getFullCorrMatrix(muAcceptErr * scaleFactor, eAcceptErr * scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(muEmErr*scaleFactor, eEmErr*scaleFactor) );//correlated
      covariance.push_back( cp.getFullCorrMatrix(muHfErr*scaleFactor, eHfErr*scaleFactor) );//correlated
      covariance.push_back( cp.getFullUncorrMatrix(muPtSmearErr*scaleFactor, ePtSmearErr*scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(muMCStatErr * scaleFactor, eMCStatErr * scaleFactor) );
      std::vector<double> combined = cp.combine(mu*scaleFactor, e*scaleFactor, covariance);

      //calculate a weighted mean
      combo[j]->SetBinContent(i, combined.at(0)/scaleFactor ); 
      combo[j]->SetBinError(i, combined.at(1)/scaleFactor );
      comboSyst[j]->SetBinContent(i, combined.at(2)/scaleFactor);
    }
    std::cout << j << std::endl;
    combo[j]->Print("All");
  }
  //**********************************************End Systematics****************************************************


  //make muon plot
  y_e->GetYaxis()->SetTitle("#frac{1}{N_{MB}} #frac{dN_{Z}}{dy}");

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  
  y_mu24->SetMarkerColor(kViolet+1);
  y_mu24->SetMarkerStyle(21);
  y_mu24->SetLineColor(kViolet+1);
  
  y_mu21->SetMarkerColor(kViolet+1);
  y_mu21->SetMarkerStyle(8);
  y_mu21->SetLineColor(kViolet+1);
  
  y_e->SetMarkerColor(kRed+1);
  y_e->SetMarkerStyle(25);
  y_e->SetLineColor(kRed+1);
  y_e->GetXaxis()->CenterTitle();
  y_e->GetYaxis()->CenterTitle();
  y_e->SetStats(0);  
  y_e->GetXaxis()->SetRangeUser(-2.4,2.4);
  y_e->GetYaxis()->SetRangeUser(0,y_e->GetMaximum()*1.4);

  y_e->Draw("p");

  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  setex2->Draw();
  EPPS16Rap->Draw("same E2");
  nCTEQ15Rap->Draw("same E2");
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
  setex1->Draw();
  y_e->Draw("same");
  //y_mu21->Draw("p same");
  y_mu24->Draw("p same");
  combo[2]->SetMarkerStyle(8);
  combo[2]->SetMarkerColor(kBlack);
  combo[2]->SetLineColor(kBlack);
  combo[2]->Draw("p same");

  TBox * eBoxy[40];
  //TBox * mu21Boxy[40];
  TBox * mu24Boxy[40];
  TBox * netBoxy[40];
  
  for(int i = 1; i<combo[2]->GetSize()-1; i++){
      helper.drawBoxAbsolute(combo[2], i , netBoxy[i], comboSyst[2]->GetBinContent(i),0.1,(Color_t)kBlack); 
      helper.drawBoxAbsolute(y_e, i , eBoxy[i], y_e->GetBinContent(i) * totalError_0_90[0][2]->GetBinContent(i) ,0.1,(Color_t)kRed+1); 
      //helper.drawBoxAbsolute(y_mu21, i , mu21Boxy[i], y_mu21->GetBinContent(i) * totalError_0_90[1][2]->GetBinContent(i),0.1,(Color_t)kViolet+1); 
  }
  for(int i = 1; i<y_mu24->GetSize()-1; i++){
      helper.drawBoxAbsolute(y_mu24, i , mu24Boxy[i], y_mu24->GetBinContent(i) * totalError_0_90[2][2]->GetBinContent(i),0.1,(Color_t)kViolet+1); 
  }

 
  y_e->Draw("p same");
  if(!doAccept) y_mu21->Draw("p same");
  y_mu24->Draw("p same");
  combo[2]->Draw("p same");

  TLegend *ly = new TLegend(0.625,0.725,0.875,0.875);
  ly->SetBorderSize(0);
  ly->SetFillStyle(0);
  ly->AddEntry((TObject*)0,"PbPb (0-90%)","");
  if(!doAccept){
    ly->AddEntry(y_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
    //ly->AddEntry(y_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
    ly->AddEntry(y_e,"Z #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
    ly->AddEntry(combo[2],"Combined (|#eta_{l}|<2.1)","p");
  }else{
    ly->AddEntry(y_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.4)","p");
    //ly->AddEntry(y_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.1)","p");
    ly->AddEntry(y_e,"Z #rightarrow e^{+}e^{-} (|y_{Z}|<2.1)","p");
    ly->AddEntry(combo[2],"Combined (|y_{Z}|<2.1)","p");
    ly->AddEntry(EPPS16Rap,"aMC@NLO + EPPS16","F");
    ly->AddEntry(nCTEQ15Rap,"aMC@NLO + nCTEQ15","F");
  }
  ly->Draw("same");
  
  c1->RedrawAxis();
  CMS_lumi(c1,0,10);
 
  if(!doAccept){ 
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty.png"); 
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty.pdf"); 
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty.C"); 
  }
  else{
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept.png"); 
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept.pdf"); 
    c1->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept.C"); 
  }
  
  //make pt plot
  pt_e->GetYaxis()->SetTitle("#frac{1}{N_{MB}} #frac{dN_{Z}}{dp_{T}} (GeV^{-1})");
  
  pt_mu24->SetMarkerColor(kBlue);
  pt_mu24->SetMarkerStyle(21);
  pt_mu24->SetLineColor(kBlue);
  
  pt_mu21->SetMarkerColor(kViolet+1);
  pt_mu21->SetMarkerStyle(8);
  pt_mu21->SetLineColor(kViolet+1);
  
  pt_e->SetMarkerColor(kRed+1);
  pt_e->SetMarkerStyle(25);
  pt_e->SetLineColor(kRed+1);
    
  TH1D * dummy = new TH1D("dummy",";p_{T} (GeV);#frac{1}{N_{MB}} #frac{dN_{Z}}{dp_{T}} (GeV^{-1})",2,0.1,200);
  dummy->SetBinContent(1,pt_mu24->GetMaximum());
  dummy->SetBinContent(2,pt_e->GetMinimum());
  dummy->SetLineColor(kWhite);
  dummy->SetMarkerColor(kWhite);
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  dummy->SetStats(0);
  dummy->Draw();
 
  setex2->Draw();
  EPPS16Pt->Draw("same E2"); 
  nCTEQ15Pt->Draw("same E2"); 
  setex1->Draw();
  pt_e->Draw("p same");
  pt_mu21->Draw("p same");
  //pt_mu24->Draw("p same");
  combo[1]->SetMarkerStyle(8);
  combo[1]->SetMarkerColor(kBlack);
  combo[1]->SetLineColor(kBlack);
  combo[1]->Draw("p same");
  
  TBox * eBoxpt[40];
  TBox * mu21Boxpt[40];
  //TBox * mu24Boxpt[40];
  TBox * netBoxpt[40];
  
  for(int i = 1; i<combo[1]->GetSize()-1; i++){
      float width = (comboSyst[1]->GetXaxis()->GetBinUpEdge(i)-comboSyst[1]->GetXaxis()->GetBinLowEdge(i)) * 0.8 / 2.0;
      helper.drawBoxAbsolute(combo[1], i , netBoxpt[i], comboSyst[1]->GetBinContent(i),width,(Color_t)kBlack); 
      helper.drawBoxAbsolute(pt_e, i , eBoxpt[i], pt_e->GetBinContent(i) * totalError_0_90[0][1]->GetBinContent(i) ,width ,(Color_t)kRed+1); 
      helper.drawBoxAbsolute(pt_mu21, i , mu21Boxpt[i], pt_mu21->GetBinContent(i) * totalError_0_90[1][1]->GetBinContent(i),width,(Color_t)kViolet+1); 
   //   helper.drawBoxAbsolute(pt_mu24, i , mu24Boxpt[i], pt_mu24->GetBinContent(i) * totalError_0_90[2][1]->GetBinContent(i),width,(Color_t)kBlue); 
  }
 

 
  pt_e->Draw("p same");
  pt_mu21->Draw("p same");
  //pt_mu24->Draw("p same");
  combo[1]->Draw("p same");

  c1->SetLogy();
  c1->SetLogx();

  TLegend *lpt = new TLegend(0.325,0.325,0.575,0.475);
  lpt->SetBorderSize(0);
  lpt->SetFillStyle(0);
  lpt->AddEntry((TObject*)0,"PbPb (0-90%)","");
  if(!doAccept){
    //lpt->AddEntry(pt_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
    lpt->AddEntry(pt_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
    lpt->AddEntry(pt_e,"Z #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
    lpt->AddEntry(combo[1],"Combined (|#eta_{l}|<2.1)","p");
  }else{
    //lpt->AddEntry(pt_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.4)","p");
    lpt->AddEntry(pt_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.1)","p");
    lpt->AddEntry(pt_e,"Z #rightarrow e^{+}e^{-} (|y_{Z}|<2.1)","p");
    lpt->AddEntry(combo[1],"Combined (|y_{Z}|<2.1)","p");
    lpt->AddEntry(EPPS16Pt,"aMC@NLO + EPPS16","F");
    lpt->AddEntry(nCTEQ15Pt,"aMC@NLO + nCTEQ15","F");
  }
  lpt->Draw("same");
  

  c1->RedrawAxis();
  CMS_lumi(c1,0,10);

  if(!doAccept){
    c1->SaveAs("plots/prettyPlots/pt_Pretty.png"); 
    c1->SaveAs("plots/prettyPlots/pt_Pretty.pdf"); 
    c1->SaveAs("plots/prettyPlots/pt_Pretty.C"); 
  }
  else{
    c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept.png"); 
    c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept.pdf"); 
    c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept.C"); 
  }

  TCanvas * cRat = new TCanvas("cRat","cRat",800,800);
  dummy->GetYaxis()->SetRangeUser(0,2.0);
  dummy->Draw(); 
  dummy->GetYaxis()->SetTitle("muon/electron channel");
  pt_mu21->Divide(pt_e);
  pt_mu21->Draw("same");
  if(!doAccept){
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty.png"); 
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty.pdf"); 
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty.C"); 
  }
  else{
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty_withAccept.png"); 
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty_withAccept.pdf"); 
    cRat->SaveAs("plots/prettyPlots/ptChannelRatio_Pretty_withAccept.C"); 
  }
  

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 7)
  {
    std::cout << "Usage: massPeakPlots <Z2EE file> <Z2mumu21 file> <Z2mumu24 file> <Z2EE syst> <Z2mumu21 syst> <Z2mumu24 syst>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu21 = argv[2];
  std::string Zmumu24 = argv[3];
  std::string ZeeSyst = argv[4];
  std::string Zmumu21Syst = argv[5];
  std::string Zmumu24Syst = argv[6];
   
  prettyPlots(Zee, Zmumu21, Zmumu24, ZeeSyst, Zmumu21Syst, Zmumu24Syst, false);
  prettyPlots(Zee, Zmumu21, Zmumu24, ZeeSyst, Zmumu21Syst, Zmumu24Syst, true);
  return 0; 
}

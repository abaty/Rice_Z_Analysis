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

  float unitScale  = 1.0;
  float unitScaleTh  = 1000.0;

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
  int offset = 1;
  for(int i = EPPS16Rap->GetSize()/2; i<EPPS16Rap->GetSize(); i++){
    EPPS16Rap->SetBinContent(i, ( EPPS16Rap->GetBinContent(i) + EPPS16Rap->GetBinContent(i-offset) )/2.0 );
    EPPS16Rap->SetBinError(i, ( EPPS16Rap->GetBinError(i) + EPPS16Rap->GetBinError(i-offset) )/2.0 );
    offset+=2;
  }

  TH1D * nCTEQ15Pt = (TH1D*)theory->Get("theorynCTEQ15_pt_Band");
  TH1D * nCTEQ15Rap = (TH1D*)theory->Get("theorynCTEQ15_rap_Band");
  offset = 1;
  for(int i = nCTEQ15Rap->GetSize()/2; i<nCTEQ15Rap->GetSize(); i++){
    nCTEQ15Rap->SetBinContent(i, ( nCTEQ15Rap->GetBinContent(i) + nCTEQ15Rap->GetBinContent(i-offset) )/2.0 );
    nCTEQ15Rap->SetBinError(i, ( nCTEQ15Rap->GetBinError(i) + nCTEQ15Rap->GetBinError(i-offset) )/2.0 );
    offset+=2;
  }
  
  TH1D * CT14Rap = (TH1D*)theory->Get("theoryCT14_rap_Band");
  TH1D * CT14Pt = (TH1D*)theory->Get("theoryCT14_pt_Band");
  offset = 1;
  for(int i = CT14Rap->GetSize()/2; i<CT14Rap->GetSize(); i++){
    CT14Rap->SetBinContent(i, ( CT14Rap->GetBinContent(i) + CT14Rap->GetBinContent(i-offset) )/2.0 );
    CT14Rap->SetBinError(i, ( CT14Rap->GetBinError(i) + CT14Rap->GetBinError(i-offset) )/2.0 );
    offset+=2;
  }

  EPPS16Pt->SetFillColor(kBlue);
  EPPS16Pt->SetLineColor(kBlue);
  EPPS16Pt->SetMarkerColorAlpha(kBlue,0);
  EPPS16Pt->SetFillStyle(1001);
  EPPS16Pt->SetFillColorAlpha(kBlue,0.5);
  EPPS16Rap->SetFillColor(kBlue);
  EPPS16Rap->SetLineColor(kBlue);
  EPPS16Rap->SetMarkerColorAlpha(kBlue,0);
  EPPS16Rap->SetFillStyle(1001);
  EPPS16Rap->SetFillColorAlpha(kBlue,0.5);
  nCTEQ15Pt->SetFillColor(kRed+1);
  nCTEQ15Pt->SetFillStyle(3354);
  nCTEQ15Pt->SetLineColor(kRed+1);
  nCTEQ15Pt->SetMarkerColorAlpha(kRed+1,0);
  nCTEQ15Rap->SetFillColor(kRed+1);
  nCTEQ15Rap->SetFillStyle(3354);
  nCTEQ15Rap->SetLineColor(kRed+1);
  nCTEQ15Rap->SetMarkerSize(0);
  nCTEQ15Rap->SetMarkerColorAlpha(kRed+1,0);
  CT14Rap->SetMarkerSize(0);
  CT14Rap->SetMarkerColorAlpha(kGreen+2,0);
  CT14Rap->SetFillColor(kGreen+2);
  CT14Rap->SetLineColor(kGreen+2);
  CT14Rap->SetFillStyle(3345); 
  CT14Pt->SetMarkerColorAlpha(kGreen+2,0);
  CT14Pt->SetFillColor(kGreen+2);
  CT14Pt->SetLineColor(kGreen+2);
  CT14Pt->SetFillStyle(3345); 
 
  //gStyle->SetErrorX(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //mumu channel first
  TFile * mu24 = TFile::Open(Zmumu24.c_str(),"read");
  TH1D * y_mu24 = (TH1D*)mu24->Get("yOS_minusAll_0_100");
  TH1D * pt_mu24 = (TH1D*)mu24->Get("unfoldedDist");

  //pt_mu24->Divide( ptCorr.correction[0] );
 
  TFile * mu21 = TFile::Open(Zmumu21.c_str(),"read");
  TH1D * y_mu21 = (TH1D*)mu21->Get("yOS_minusAll_0_100");
  TH1D * pt_mu21 = (TH1D*)mu21->Get("unfoldedDist");
  
  //pt_mu21->Divide( ptCorr.correction[2] );
  
  TFile * e = TFile::Open(Zee.c_str(),"read");
  TH1D * y_e = (TH1D*)e->Get("yOS_minusAll_0_100");
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
  combo[1] = (TH1D*) pt_e->Clone("combined_pT_0_100"); 
  combo[2] = (TH1D*) y_e->Clone("combined_y_0_100"); 
  TH1D * comboSyst[3];
  comboSyst[1] = (TH1D*) pt_e->Clone("combinedSyst_pT_0_100"); 
  comboSyst[2] = (TH1D*) y_e->Clone("combinedSyst_y_0_100"); 

  for(int i = 0; i<3; i++){
    for(int j = 1; j<3; j++){
      std::cout << i << " " << j << std::endl;
      efficiencyError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_efficiencyError_0_100",helper.name.at(j).c_str()));
      emError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_emError_0_100",helper.name.at(j).c_str()));
      hfError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_hfError_0_100",helper.name.at(j).c_str()));
      if(j==1) ptSmearError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_ptSmearError_0_100",helper.name.at(j).c_str()));
      mcStatError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_mcStatError_0_100",helper.name.at(j).c_str()));
      acceptError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_acceptError_0_100",helper.name.at(j).c_str()));
      totalError_0_90[i][j] = (TH1D*) systFile[i]->Get(Form("%s_totalError_0_100",helper.name.at(j).c_str()));
    }
  }
  //plots here

  
  //scaling
  y_mu24->Scale(unitScale/(s.muLumi));
  y_mu21->Scale(unitScale/(s.muLumi));
  y_e->Scale(unitScale/(s.eLumi));
  pt_mu24->Scale(1.0/(s.muLumi));
  pt_mu21->Scale(1.0/(s.muLumi));
  pt_e->Scale(1.0/(s.eLumi));

  if(doAccept){
    pt_e->Divide(acceptE[0]);
    pt_mu21->Divide(acceptMu21[0]);
    pt_mu24->Divide(acceptMu24[0]);
    y_e->Divide(acceptE[1]);
    y_mu21->Divide(acceptMu21[1]);
    y_mu24->Divide(acceptMu24[1]);
  }

  y_e->Print("All");

  //combining positive and negative sides
  offset = 1;
  std::cout << "starting combining plus and minus sides" << std::endl;
  for(int i = y_e->GetSize()/2; i<y_e->GetSize(); i++){
    float sig1 = y_e->GetBinError(i);
    float sig2 = y_e->GetBinError(i-offset);
    if(sig1==0 || sig2==0) continue;
    float totalWeight = 1.0/(sig1*sig1) + 1.0/(sig2*sig2);
    //using BLUE for stat errors (systematics are all correlated except MC stats which are very small)
    y_e->SetBinContent(i, ( y_e->GetBinContent(i)/(sig1*sig1) + y_e->GetBinContent(i-offset)/(sig2*sig2) )/ totalWeight);
    y_e->SetBinError(i, 1.0/TMath::Sqrt(totalWeight) );
    offset+=2;
  }
  std::cout << "done starting combining plus and minus sides for electrons" << std::endl;
  y_e->Print("All");
  offset = 1;
  for(int i = y_mu21->GetSize()/2; i<y_mu21->GetSize(); i++){
    float sig1 = y_mu21->GetBinError(i);
    float sig2 = y_mu21->GetBinError(i-offset);
    float totalWeight = 1.0/(sig1*sig1) + 1.0/(sig2*sig2);
    //using BLUE for stat errors (systematics are all correlated except MC stats which are very small)
    y_mu21->SetBinContent(i, ( y_mu21->GetBinContent(i)/(sig1*sig1) + y_mu21->GetBinContent(i-offset)/(sig2*sig2) )/ totalWeight);
    y_mu21->SetBinError(i, 1.0/TMath::Sqrt(totalWeight) );
    offset+=2;
  }
  offset = 1;
  for(int i = y_mu24->GetSize()/2; i<y_mu24->GetSize(); i++){
    float sig1 = y_mu24->GetBinError(i);
    float sig2 = y_mu24->GetBinError(i-offset);
    float totalWeight = 1.0/(sig1*sig1) + 1.0/(sig2*sig2);
    //using BLUE for stat errors (systematics are all correlated except MC stats which are very small)
    y_mu24->SetBinContent(i, ( y_mu24->GetBinContent(i)/(sig1*sig1) + y_mu24->GetBinContent(i-offset)/(sig2*sig2) )/ totalWeight);
    y_mu24->SetBinError(i, 1.0/TMath::Sqrt(totalWeight) );
    offset+=2;
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
  y_e->GetYaxis()->SetTitle("#frac{d#sigma_{Z}}{d|y_{Z}|} (#mub)");

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  
  y_mu24->SetMarkerColor(kBlue);
  y_mu24->SetMarkerStyle(21);
  y_mu24->SetLineColor(kBlue);
  
  y_mu21->SetMarkerColor(kBlue);
  y_mu21->SetMarkerStyle(8);
  y_mu21->SetLineColor(kBlue);
  
  y_e->SetMarkerColor(kRed+1);
  y_e->SetMarkerStyle(25);
  y_e->SetLineColor(kRed+1);
  y_e->GetXaxis()->CenterTitle();
  y_e->GetYaxis()->CenterTitle();
  y_e->SetStats(0);  
  y_e->GetXaxis()->SetRangeUser(-2.4,2.4);
  y_e->GetYaxis()->SetRangeUser(0,y_e->GetMaximum()*1.4);

  y_e->GetYaxis()->SetTitleSize(0.05);
  y_e->GetYaxis()->SetTitleOffset(1.6);
  y_e->GetYaxis()->SetLabelSize(0.05);
  y_e->GetYaxis()->SetLabelOffset(0.005);
  y_e->GetXaxis()->SetTitleSize(0.06);
  y_e->GetXaxis()->SetTitleOffset(1);
  y_e->GetXaxis()->SetLabelSize(0.04);
  y_e->GetXaxis()->SetLabelOffset(0.005);

  y_e->Draw("p");

  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  setex2->Draw();
  //float TAA090 = 6.274;
  //float TAA0100 = 5.649;
  //EPPS16Rap->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  //nCTEQ15Rap->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  //CT14Rap->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  EPPS16Rap->Scale(unitScaleTh);//need to convert to ub
  nCTEQ15Rap->Scale(unitScaleTh);//need to convert to ub
  CT14Rap->Scale(unitScaleTh);//need to convert to ub
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
      //helper.drawBoxAbsolute(y_mu21, i , mu21Boxy[i], y_mu21->GetBinContent(i) * totalError_0_90[1][2]->GetBinContent(i),0.1,(Color_t)kBlue); 
  }
  for(int i = 1; i<y_mu24->GetSize()-1; i++){
      helper.drawBoxAbsolute(y_mu24, i , mu24Boxy[i], y_mu24->GetBinContent(i) * totalError_0_90[2][2]->GetBinContent(i),0.1,(Color_t)kBlue); 
  }

 
  y_e->Draw("p same");
  if(!doAccept) y_mu21->Draw("p same");
  y_mu24->Draw("p same");
  combo[2]->Draw("p same");

  TLegend *ly = new TLegend(0.225,0.225,0.705,0.5);
  ly->SetBorderSize(0);
  ly->SetFillStyle(0);
  ly->AddEntry((TObject*)0,"0-100%","");
  if(!doAccept){
    ly->AddEntry(y_mu24,"Z/#gamma* #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
    //ly->AddEntry(y_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
    ly->AddEntry(y_e,"Z/#gamma* #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
    ly->AddEntry(combo[2],"Combined (|#eta_{l}|<2.1)","p");
  }else{
    ly->AddEntry(y_mu24,"Z/#gamma* #rightarrow #mu^{+}#mu^{-}","p");
    //ly->AddEntry(y_mu21,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.1)","p");
    ly->AddEntry(y_e,"Z/#gamma* #rightarrow e^{+}e^{-}","p");
    ly->AddEntry(combo[2],"Combined","p");
    ly->AddEntry(EPPS16Rap,"MG5_aMC@NLO + CT14 + EPPS16","F");
    ly->AddEntry(nCTEQ15Rap,"MG5_aMC@NLO + nCTEQ15","F");
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
 
  TCanvas * c2 = new TCanvas("c2","c2",800,500);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.025);
  c2->SetBottomMargin(0.16);  
  c2->SetTopMargin(0.1);  
  y_e->SetLineColor(kWhite);
  y_e->SetMarkerColor(kWhite);
  y_e->GetYaxis()->SetRangeUser(2.2,5.7);
  y_e->GetYaxis()->SetTitleSize(0.075);
  y_e->GetYaxis()->SetTitleOffset(0.75);
  y_e->GetYaxis()->SetLabelSize(0.07);
  y_e->GetXaxis()->SetTitleSize(0.08);
  y_e->GetXaxis()->SetTitleOffset(0.85);
  y_e->GetXaxis()->SetLabelSize(0.07);
  y_e->GetYaxis()->SetNdivisions(5,5,0,true);
  y_e->GetXaxis()->SetTitle("|y_{Z}|");
  y_e->GetXaxis()->SetNdivisions(408,false);
  y_e->GetXaxis()->SetRangeUser(0,2.5);
  y_e->Draw("p");

  for(int i = 0; i<y_mu24->GetSize(); i++){
    if(i!=1 && (i != y_mu24->GetSize()-2)){
      y_mu24->SetBinContent(i,0);
      y_mu24->SetBinError(i,0);
    }
  }

  y_mu24->SetMarkerStyle(25);
  y_mu24->SetMarkerColor(kBlack);
  y_mu24->SetLineColor(kBlack);
  setex2->Draw();
  EPPS16Rap->Draw("same E2");
  nCTEQ15Rap->Draw("same E2");
  CT14Rap->Draw("same E2");
  combo[2]->Draw("p same");
  y_mu24->Draw("p same");
  for(int i = combo[2]->GetSize()/2; i<combo[2]->GetSize()-1; i++){
      helper.drawBoxAbsolute(combo[2], i , netBoxy[i], comboSyst[2]->GetBinContent(i),0.145,(Color_t)kBlack); 
  }
  for(int i = 1; i<y_mu24->GetSize()-1; i++){
      if(i==1 || i==y_mu24->GetSize()-2){
        helper.drawBoxAbsolute(y_mu24, i , mu24Boxy[i], y_mu24->GetBinContent(i) * totalError_0_90[2][2]->GetBinContent(i),0.145,(Color_t)kBlack); 
      }
  }
  TLegend *ly2 = new TLegend(0.165,0.19,0.765,0.565);
  ly2->SetBorderSize(0);
  ly2->SetFillStyle(0);
  ly2->SetTextSize(0.065);
  ly2->AddEntry(combo[2],"Z/#gamma* #rightarrow l^{+}l^{-}","ple");
  ly2->AddEntry(y_mu24,"Z/#gamma* #rightarrow #mu^{+}#mu^{-}","ple");
  ly2->AddEntry(CT14Rap,"MG5_aMC@NLO + CT14","F");
  ly2->AddEntry(nCTEQ15Rap,"MG5_aMC@NLO + nCTEQ15","F");
  ly2->AddEntry(EPPS16Rap,"MG5_aMC@NLO + CT14 + EPPS16","F");
  ly2->Draw("same");
  TLegend * ly3 = new TLegend(0.58,0.8,0.88,0.875);
  ly3->SetTextSize(0.065);
  ly3->SetBorderSize(0);
  ly3->SetFillStyle(0);
  ly3->AddEntry((TObject*)0,"60 < m_{ll} < 120 GeV","");
  ly3->Draw("same");
  
  c2->RedrawAxis();
  CMS_lumi(c2,0,0,1.5, true, true, false, true);

  std::cout << "Rapidity Compatibility" << std::endl;
  float maxDeviationY = -99;
  for(int i = 1; i<combo[2]->GetSize()-1; i++){
    float e = y_e->GetBinContent(i);
    float mu = y_mu21->GetBinContent(i);
    float difference = TMath::Abs(e-mu);
 
    float eStat = y_e->GetBinError(i);
    float muStat = y_mu21->GetBinError(i);
    float eSyst = y_e->GetBinContent(i) * totalError_0_90[0][2]->GetBinContent(i);
    float muSyst = y_mu21->GetBinContent(i) * totalError_0_90[1][2]->GetBinContent(i);
    float err = TMath::Sqrt( eStat * eStat + muStat * muStat + eSyst * eSyst + muSyst * muSyst );
    std::cout << i << " " << difference/err << " sigma" << std::endl;

    if(difference/err > maxDeviationY) maxDeviationY = difference/err;
  }
  std::cout << "Max deviation: " << maxDeviationY << std::endl;

  c2->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept_PRL.png"); 
  c2->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept_PRL.pdf"); 
  c2->SaveAs("plots/prettyPlots/rapidity_Pretty_withAccept_PRL.C"); 


  delete c1;
  c1 = new TCanvas("c1","c1",800,800);
  TPad * p1 = new TPad("p1","p1",0,0.2,1,1,0);
  TPad * p2 = new TPad("p2","p2",0,0,1,0.2,0);
  c1->SetLineWidth(0);
  p1->SetBottomMargin(0);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.05);
  p1->SetTopMargin(0.07);
  p1->SetBorderSize(0);
  p1->Draw();
  p2->SetTopMargin(0);
  p2->SetLeftMargin(0.15);
  p2->SetRightMargin(0.05);
  p2->SetBottomMargin(0.3);
  p2->SetBorderSize(0);
  p2->Draw();
  p1->cd();
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
 
  //make pt plot
  pt_e->GetYaxis()->SetTitle("#frac{d#sigma_{Z}}{dp_{T}^{Z}} (#mub/GeV)");
  
  pt_mu24->SetMarkerColor(kBlue);
  pt_mu24->SetMarkerStyle(21);
  pt_mu24->SetLineColor(kBlue);
  
  pt_mu21->SetMarkerColor(kBlue);
  pt_mu21->SetMarkerStyle(8);
  pt_mu21->SetLineColor(kBlue);
  
  pt_e->SetMarkerColor(kRed+1);
  pt_e->SetMarkerStyle(25);
  pt_e->SetLineColor(kRed+1);
    
  TH1D * dummy = new TH1D("dummy",";p_{T}^{Z} (GeV);#frac{d#sigma_{Z}}{dp_{T}^{Z}} (#mub/GeV)",2,0.1,200);
  dummy->SetBinContent(1,pt_mu24->GetMaximum());
  dummy->SetBinContent(2,pt_e->GetMinimum());
  dummy->SetLineColor(kWhite);
  dummy->SetMarkerColor(kWhite);
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  dummy->GetYaxis()->SetTitleSize(0.05);
  dummy->GetYaxis()->SetTitleOffset(1.27);
  dummy->GetYaxis()->SetLabelSize(0.05);
  dummy->GetYaxis()->SetLabelOffset(0.005);
  dummy->SetStats(0);
  dummy->DrawCopy();
 
  setex2->Draw();
  
  //EPPS16Pt->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  //nCTEQ15Pt->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  //CT14Pt->Scale(unitScale * s.sigmaMBPbPb * TMath::Power(10,6) * TAA0100/TAA090);//need to convert to ub
  EPPS16Pt->Scale(unitScaleTh);//need to convert to ub
  nCTEQ15Pt->Scale(unitScaleTh);//need to convert to ub
  CT14Pt->Scale(unitScaleTh);//need to convert to ub
  
  EPPS16Pt->Draw("same E2"); 
  nCTEQ15Pt->Draw("same E2"); 
  CT14Pt->Draw("same E2");
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
      float width = (comboSyst[1]->GetXaxis()->GetBinCenter(i)-comboSyst[1]->GetXaxis()->GetBinLowEdge(i));
      float width2 = (comboSyst[1]->GetXaxis()->GetBinUpEdge(i)-comboSyst[1]->GetXaxis()->GetBinCenter(i));
      helper.drawBoxAbsolute(combo[1], i , netBoxpt[i], comboSyst[1]->GetBinContent(i),width,(Color_t)kBlack,true,0,width2); 
      helper.drawBoxAbsolute(pt_e, i , eBoxpt[i], pt_e->GetBinContent(i) * totalError_0_90[0][1]->GetBinContent(i) ,width ,(Color_t)kRed+1, true, 0 , width2); 
      helper.drawBoxAbsolute(pt_mu21, i , mu21Boxpt[i], pt_mu21->GetBinContent(i) * totalError_0_90[1][1]->GetBinContent(i),width,(Color_t)kBlue, true, 0, width2); 
   //   helper.drawBoxAbsolute(pt_mu24, i , mu24Boxpt[i], pt_mu24->GetBinContent(i) * totalError_0_90[2][1]->GetBinContent(i),width,(Color_t)kBlue); 
  }
 
  pt_e->Draw("p same");
  pt_mu21->Draw("p same");
  //pt_mu24->Draw("p same");
  combo[1]->Draw("p same");

  p1->SetLogy();
  p1->SetLogx();


  TLegend *lpt = new TLegend(0.225,0.125,0.675,0.475);
  lpt->SetBorderSize(0);
  lpt->SetFillStyle(0);
  lpt->AddEntry((TObject*)0,"0-100%, |y_{Z}|<2.1","");
  if(!doAccept){
    //lpt->AddEntry(pt_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.4)","p");
    lpt->AddEntry(pt_mu21,"Z/#gamma* #rightarrow #mu^{+}#mu^{-} (|#eta_{#mu}|<2.1)","p");
    lpt->AddEntry(pt_e,"Z/#gamma* #rightarrow e^{+}e^{-} (|#eta_{e}|<2.1)","p");
    lpt->AddEntry(combo[1],"Combined (|#eta_{l}|<2.1)","p");
  }else{
    //lpt->AddEntry(pt_mu24,"Z #rightarrow #mu^{+}#mu^{-} (|y_{Z}|<2.4)","p");
    lpt->AddEntry(pt_mu21,"Z/#gamma* #rightarrow #mu^{+}#mu^{-}","p");
    lpt->AddEntry(pt_e,"Z/#gamma* #rightarrow e^{+}e^{-}","p");
    lpt->AddEntry(combo[1],"Combined","p");
    lpt->AddEntry(EPPS16Pt,"MG5_aMC@NLO + CT14 + EPPS16","F");
    lpt->AddEntry(nCTEQ15Pt,"MG5_aMC@NLO + nCTEQ15","F");
  }
  lpt->Draw("same");
  

  p1->RedrawAxis();
  CMS_lumi(p1,0,10,1.8, true, true, false, true);
  
  p2->cd();
  p2->SetLogx();
  dummy->GetYaxis()->SetTitle("#frac{MC}{Data}");
  dummy->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
  dummy->GetYaxis()->SetTitleOffset(0.4);
  dummy->GetYaxis()->SetRangeUser(0.4,1.6);
  dummy->GetYaxis()->SetTitleSize(0.21);
  dummy->GetYaxis()->SetTitleOffset(0.31);
  dummy->GetXaxis()->SetTitleSize(0.18);
  dummy->GetXaxis()->SetTitleOffset(0.7);
  dummy->GetYaxis()->SetLabelSize(0.19);
  dummy->GetYaxis()->SetLabelOffset(0.005);
  dummy->GetXaxis()->SetLabelSize(0.21);
  dummy->GetXaxis()->SetLabelOffset(-0.053);
  dummy->Draw();


  TH1D * datErr = (TH1D*)EPPS16Pt->Clone("datErr");
  for(int i = 1; i<datErr->GetSize()-1; i++){
    float err = TMath::Power( TMath::Power(combo[1]->GetBinError(i)/combo[1]->GetBinContent(i),2) + TMath::Power(comboSyst[1]->GetBinContent(i)/combo[1]->GetBinContent(i),2) , 0.5 );
    datErr->SetBinContent(i,1);
    datErr->SetBinError(i,err);
  }
  datErr->Print("All");
  datErr->SetFillStyle(1001);
  datErr->SetFillColor(kGray);
  datErr->SetLineColor(kGray);
  setex2->Draw();
  datErr->Draw("same E2");
  setex1->Draw();
  
  TLine * line1;
  line1 = new TLine(EPPS16Pt->GetXaxis()->GetBinLowEdge(1),1,EPPS16Pt->GetXaxis()->GetBinUpEdge(EPPS16Pt->GetSize()-2),1);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);
  line1->Draw("same");

  
  TH1D * ratioEPPS16 = (TH1D*)EPPS16Pt->Clone("ratioEPPS16");
  ratioEPPS16->Divide(combo[1]);
  ratioEPPS16->SetLineColor(kBlue);
  ratioEPPS16->SetLineWidth(2);
  ratioEPPS16->SetMarkerColor(kBlue);
  ratioEPPS16->SetMarkerStyle(8);
  ratioEPPS16->Draw("same");
  TH1D * ratioNCTEQ15 = (TH1D*)nCTEQ15Pt->Clone("ratioNCTEQ15");
  ratioNCTEQ15->Divide(combo[1]);
  ratioNCTEQ15->SetLineColor(kRed+1);
  ratioNCTEQ15->SetMarkerColor(kRed+1);
  ratioNCTEQ15->SetMarkerStyle(32);
  ratioNCTEQ15->SetMarkerSize(1.2);
  ratioNCTEQ15->Draw("same p");
  TH1D * ratioCT14 = (TH1D*)CT14Pt->Clone("ratioCT14");
  ratioCT14->Divide(combo[1]);
  ratioCT14->SetLineColor(kGreen+2);
  ratioCT14->SetMarkerColor(kGreen+2);
  ratioCT14->SetMarkerStyle(28);
  ratioCT14->SetMarkerSize(1);
  ratioCT14->Draw("same p");
  ratioCT14->Print("All"); 
  p2->RedrawAxis();
  
  TLegend *lptR = new TLegend(0.175,0.3,0.675,0.6);
  lptR->SetBorderSize(0);
  lptR->SetFillStyle(0);
  lptR->SetNColumns(2);
  lptR->AddEntry(ratioNCTEQ15,"nCTEQ15","p");
  lptR->AddEntry(ratioCT14,"CT14","p");
  lptR->AddEntry(ratioEPPS16,"CT14 + EPPS16","p");
  lptR->Draw("same");

  std::cout << "pT Compatibility" << std::endl;
  float maxDeviationPT = -99;
  for(int i = 1; i<combo[1]->GetSize()-1; i++){
    float e = pt_e->GetBinContent(i);
    float mu = pt_mu21->GetBinContent(i);
    float difference = TMath::Abs(e-mu);
 
    float eStat = pt_e->GetBinError(i);
    float muStat = pt_mu21->GetBinError(i);
    float eSyst = pt_e->GetBinContent(i) * totalError_0_90[0][1]->GetBinContent(i);
    float muSyst = pt_mu21->GetBinContent(i) * totalError_0_90[1][1]->GetBinContent(i);
    float err = TMath::Sqrt( eStat * eStat + muStat * muStat + eSyst * eSyst + muSyst * muSyst );
    std::cout << i << " " << difference/err << " sigma" << std::endl;
    
    if(difference/err > maxDeviationPT) maxDeviationPT = difference/err;
  }
  std::cout << "Max deviation: " << maxDeviationPT << std::endl;


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
  
  delete p1;
  delete p2;
  delete c1;
  c1 = new TCanvas("c1","c1",800,500);
  p1 = new TPad("p1","p1",0,0.3,1,1,0);
  p2 = new TPad("p2","p2",0,0,1,0.3,0);
  c1->SetLineWidth(0);
  p1->SetBottomMargin(0);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.001);
  p1->SetTopMargin(0.09);
  p1->SetBorderSize(0);
  p1->Draw();
  p2->SetTopMargin(0);
  p2->SetLeftMargin(0.15);
  p2->SetRightMargin(0.001);
  p2->SetBottomMargin(0.35);
  p2->SetBorderSize(0);
  p2->Draw();
  p1->cd();
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  dummy->SetBinContent(1,pt_mu24->GetMaximum());
  dummy->SetBinContent(2,pt_e->GetMinimum());
  dummy->SetLineColor(kWhite);
  dummy->SetMarkerColor(kWhite);
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  dummy->GetYaxis()->SetTitleSize(0.085);
  dummy->GetYaxis()->SetTitleOffset(0.75);
  dummy->GetYaxis()->SetLabelSize(0.065);
  dummy->GetYaxis()->SetLabelOffset(0.00);
  dummy->GetYaxis()->SetRangeUser(0.0004*pt_mu24->GetMaximum() , pt_mu24->GetMaximum()*1.4 ); 
  dummy->GetYaxis()->SetTitle("#frac{d#sigma_{Z}}{dp_{T}^{Z}} (#mub/GeV)");
  dummy->DrawCopy();
 
  setex2->Draw();
  EPPS16Pt->Draw("same E2"); 
  nCTEQ15Pt->Draw("same E2");
  CT14Pt->Draw("same E2"); 
  for(int i = 1; i<combo[1]->GetSize()-1; i++){
      float width = (comboSyst[1]->GetXaxis()->GetBinCenter(i)-comboSyst[1]->GetXaxis()->GetBinLowEdge(i));
      float width2 = (comboSyst[1]->GetXaxis()->GetBinUpEdge(i)-comboSyst[1]->GetXaxis()->GetBinCenter(i));
      helper.drawBoxAbsolute(combo[1], i , netBoxpt[i], comboSyst[1]->GetBinContent(i),width,(Color_t)kBlack,true,0,width2); 
  }
  setex1->Draw();
  combo[1]->Draw("p same");
  p1->SetLogy();
  p1->SetLogx();

  TLegend *lpt4 = new TLegend(0.65,0.7,0.85,0.88);
  lpt4->SetBorderSize(0);
  lpt4->SetFillStyle(0);
  lpt4->SetTextSize(0.08);
  lpt4->AddEntry((TObject*)0,"60 < m_{ll} < 120 GeV","");
  lpt4->AddEntry((TObject*)0,"               |y_{Z}| < 2.1","");
  lpt4->Draw("same");

  TLegend *lpt2 = new TLegend(0.175,0.04,0.775,0.665);
  lpt2->SetBorderSize(0);
  lpt2->SetFillStyle(0);
  lpt2->SetTextSize(0.09);
  lpt2->AddEntry(combo[1],"Z/#gamma* #rightarrow l^{+}l^{-}","lep");
  lpt2->AddEntry(CT14Pt,"MG5_aMC@NLO + CT14","F");
  lpt2->AddEntry(nCTEQ15Pt,"MG5_aMC@NLO + nCTEQ15","F");
  lpt2->AddEntry(EPPS16Pt,"MG5_aMC@NLO + CT14 + EPPS16","F");
  lpt2->Draw("same");
  
  p1->RedrawAxis();
  CMS_lumi(p1,0,0,2.0, true, true, false, true,true);
  
  p2->cd();
  p2->SetLogx();
  dummy->GetYaxis()->SetTitle("#frac{MC}{Data}");
  dummy->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
  dummy->GetYaxis()->SetTitleOffset(0.4);
  dummy->GetYaxis()->SetRangeUser(0.15,1.65);
  dummy->GetYaxis()->SetTitleSize(0.22);
  dummy->GetYaxis()->SetTitleOffset(0.3);
  dummy->GetXaxis()->SetTitleSize(0.18);
  dummy->GetXaxis()->SetTitleOffset(0.8);
  dummy->GetYaxis()->SetLabelSize(0.16);
  dummy->GetYaxis()->SetLabelOffset(0.005);
  dummy->GetXaxis()->SetLabelSize(0.17);
  dummy->GetXaxis()->SetLabelOffset(0.0);
  dummy->Draw();
  
  setex2->Draw();
  datErr->Draw("same E2");
  setex1->Draw();
  line1->Draw("same");

  ratioEPPS16->Draw("same");
  ratioNCTEQ15->Draw("same p");
  ratioCT14->Draw("same p");
  p2->RedrawAxis();
  
  TLegend *lptR2 = new TLegend(0.16,0.35,0.99,0.6);
  lptR2->SetBorderSize(0);
  lptR2->SetFillStyle(0);
  lptR2->SetNColumns(4);
  lptR2->SetTextSize(0.17);
  lptR2->AddEntry(ratioCT14,"CT14","p");
  lptR2->AddEntry(ratioEPPS16,"CT14 + EPPS16","p");
  lptR2->AddEntry(ratioNCTEQ15,"nCTEQ15","p");
  lptR2->AddEntry(datErr,"Syst. uncertainty","f");
  lptR2->Draw("same");

  c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept_PRL.png"); 
  c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept_PRL.pdf"); 
  c1->SaveAs("plots/prettyPlots/pt_Pretty_withAccept_PRL.C"); 

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

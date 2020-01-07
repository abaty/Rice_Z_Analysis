#include "TExec.h"
#include "include/CMS_lumi.C"
#include "include/centralityTool.h"
#include "include/combinePoints.h"
#include "include/HistNameHelper.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void plotMassPeaks(std::string Zee, std::string Zmumu21, std::string Zmumu24, std::string ZeeSyst, std::string Zmumu21Syst, std::string Zmumu24Syst, bool doAccept){
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  CombinePoints cp = CombinePoints();
  HistNameHelper helper = HistNameHelper();

  TFile * HGPythia = TFile::Open("resources/HGPythia_hfSum.root","read");
  TH1D * hgPythia = (TH1D*)HGPythia->Get("RAA_hf_rebin");
 
  TH1D * acceptE[1];
  TH1D * acceptMu21[1];
  TH1D * acceptMu24[1];
  if(doAccept){
    TFile * f_acceptE = TFile::Open("resources/Z2ee_EfficiencyMC_0.root");
    acceptE[0] = (TH1D*) f_acceptE->Get("accept21_yields_ratio_0");
    TFile * f_acceptMu = TFile::Open("resources/Z2mumu_Efficiencies.root");
    acceptMu21[0] = (TH1D*) f_acceptMu->Get("accept21_yields_ratio_0");
    acceptMu24[0] = (TH1D*) f_acceptMu->Get("accept24_yields_ratio_0");

    for(int i = 0; i<1; i++){
      for(int j = 0; j<acceptE[i]->GetSize(); j++)  acceptE[i]->SetBinError(j,0);
      for(int j = 0; j<acceptMu21[i]->GetSize(); j++)  acceptMu21[i]->SetBinError(j,0);
      for(int j = 0; j<acceptMu24[i]->GetSize(); j++)  acceptMu24[i]->SetBinError(j,0);
    }
  }

  gStyle->SetErrorX(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TH1D * mu24[nBins];
  TH1D * mu21[nBins];
  TH1D * e[nBins];
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZmumuFile21 = TFile::Open(Zmumu21.c_str(),"read");
  TFile * ZmumuFile24 = TFile::Open(Zmumu24.c_str(),"read");
   
  TFile * fZeeSyst = TFile::Open(ZeeSyst.c_str(),"read"); 
  TFile * fZmumu21Syst = TFile::Open(Zmumu21Syst.c_str(),"read"); 
  TFile * fZmumu24Syst = TFile::Open(Zmumu24Syst.c_str(),"read"); 

  TH1D * efficiencyError[nBins][3];
  TH1D * emError[nBins][3];	
  TH1D * hfError[nBins][3];	
  TH1D * acceptError[nBins][3];
  TH1D * totalError[nBins][3];	
  
  TH1D * combo[nBins];
  TH1D * comboSyst[nBins];

  for(int i = 0; i<nBins; i++){
    mu24[i] = (TH1D*)ZmumuFile24->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    mu21[i] = (TH1D*)ZmumuFile21->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    e[i] = (TH1D*)ZeeFile->Get(Form("yieldOS_minusAll_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    combo[i] = (TH1D*) e[i]->Clone(Form("yieldCombo_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i))); 
    comboSyst[i] = (TH1D*) e[i]->Clone(Form("yieldComboSyst_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i))); 
  
    efficiencyError[i][0] = (TH1D*)fZeeSyst->Get(Form("yield_efficiencyError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    efficiencyError[i][1] = (TH1D*)fZmumu21Syst->Get(Form("yield_efficiencyError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    efficiencyError[i][2] = (TH1D*)fZmumu24Syst->Get(Form("yield_efficiencyError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    acceptError[i][0] = (TH1D*)fZeeSyst->Get(Form("yield_acceptError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    acceptError[i][1] = (TH1D*)fZmumu21Syst->Get(Form("yield_acceptError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    acceptError[i][2] = (TH1D*)fZmumu24Syst->Get(Form("yield_acceptError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    emError[i][0] = (TH1D*)fZeeSyst->Get(Form("yield_emError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    emError[i][1] = (TH1D*)fZmumu21Syst->Get(Form("yield_emError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    emError[i][2] = (TH1D*)fZmumu24Syst->Get(Form("yield_emError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    hfError[i][0] = (TH1D*)fZeeSyst->Get(Form("yield_hfError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    hfError[i][1] = (TH1D*)fZmumu21Syst->Get(Form("yield_hfError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    hfError[i][2] = (TH1D*)fZmumu24Syst->Get(Form("yield_hfError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    
    totalError[i][0] = (TH1D*)fZeeSyst->Get(Form("yield_totalError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    totalError[i][1] = (TH1D*)fZmumu21Syst->Get(Form("yield_totalError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    totalError[i][2] = (TH1D*)fZmumu24Syst->Get(Form("yield_totalError_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
  
    //scale by lumi and nMB
    mu24[i]->Scale(s.netLumi/(s.muLumi * s.Nmb));
    mu21[i]->Scale(s.netLumi/(s.muLumi * s.Nmb));
    e[i]->Scale(s.netLumi/(s.eLumi * s.Nmb));
  }

  if(doAccept){
    for(int i = 0; i<nBins; i++){
      e[i]->Divide(acceptE[0]);
      mu21[i]->Divide(acceptMu21[0]);
      mu24[i]->Divide(acceptMu24[0]);
    }
  }

  //combination
  for(int j = 0; j<nBins; j++){
    TH1D * tempHistE, *tempHistMu;
    tempHistE = e[j];
    tempHistMu = mu21[j];

    tempHistE->Print("All");
    tempHistMu->Print("All");

    for(int i = 0; i<combo[j]->GetSize()+2; i++){
      float e = tempHistE->GetBinContent(i);
      float eStatErr = tempHistE->GetBinError(i);
      float eEffErr = efficiencyError[j][0]->GetBinContent(i) * e;
      float eEmErr = emError[j][0]->GetBinContent(i) * e;
      float eHfErr = hfError[j][0]->GetBinContent(i) * e;
      float eChargeSwapErr = 0.005 * e;
      float eAcceptErr = acceptError[j][0]->GetBinContent(i) * e;
      
      float mu = tempHistMu->GetBinContent(i);
      float muStatErr = tempHistMu->GetBinError(i);
      float muEffErr = efficiencyError[j][1]->GetBinContent(i) * mu;
      float muEmErr = emError[j][1]->GetBinContent(i) * mu;
      float muHfErr = hfError[j][1]->GetBinContent(i) * mu;
      float muAcceptErr = acceptError[j][1]->GetBinContent(i) * mu;

      double scaleFactor = 10000000;     
 
      std::vector< TMatrixD > covariance;
      covariance.push_back( cp.getFullUncorrMatrix(muStatErr*scaleFactor, eStatErr*scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(0, eChargeSwapErr*scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(muEffErr*scaleFactor, eEffErr*scaleFactor) );
      covariance.push_back( cp.getFullCorrMatrix(muAcceptErr*scaleFactor, eAcceptErr * scaleFactor) );
      covariance.push_back( cp.getFullUncorrMatrix(muEmErr*scaleFactor, eEmErr*scaleFactor) );//correlated
      covariance.push_back( cp.getFullCorrMatrix(muHfErr*scaleFactor, eHfErr*scaleFactor) );//correlated
      std::vector<double> combined = cp.combine(mu*scaleFactor, e*scaleFactor, covariance);

      //calculate a weighted mean
      combo[j]->SetBinContent(i, combined.at(0)/scaleFactor ); 
      combo[j]->SetBinError(i, combined.at(1)/scaleFactor );
      comboSyst[j]->SetBinContent(i, combined.at(2)/scaleFactor);
    }
  }



  const char * labels[10] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%", "50-70%", "70-90%", "0-90%","pp"};
  float TAA[9] = {25.82, 20.36, 14.36, 8.778, 5.099, 2.778, 1.010, 0.1674, 6.274};
  float TAARelErr[9] = {0.018, 0.021, 0.023, 0.027, 0.033, 0.041, 0.045, 0.046, 0.022};
  float scaleFactor_mumu[9];
  for(int i = 0; i<9; i++){
    scaleFactor_mumu[i] = 1.0/TAA[i];
    if(i ==0 || i==1) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_mumu[i] = scaleFactor_mumu[i] * 4.5;
  }
  float scaleFactor_ee[9];
  for(int i = 0; i<9; i++){
    scaleFactor_ee[i] = 1.0/TAA[i];
    if(i ==0 || i==1) scaleFactor_ee[i] = scaleFactor_ee[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_ee[i] = scaleFactor_ee[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_ee[i] = scaleFactor_ee[i] * 4.5;
  }
  
  float scaleFactor_Combo[9];
  for(int i = 0; i<9; i++){
    scaleFactor_Combo[i] = 1.0/TAA[i];
    if(i ==0 || i==1) scaleFactor_Combo[i] = scaleFactor_Combo[i] * 18.0;
    if(i == 2 || i==3 || i==4 || i==5 ) scaleFactor_Combo[i] = scaleFactor_Combo[i] * 9.0;
    if(i == 6 || i==7) scaleFactor_Combo[i] = scaleFactor_Combo[i] * 4.5;
  }

 // gStyle->SetErrorX(0);

  int binMap[9] = {0,1,2,3,4,5,15,16,25};

  TH1D * yieldPlot_mumu24 = new TH1D("yieldPlot_mumu24","",10,0,10);
  yieldPlot_mumu24->SetBinContent(1,mu24[binMap[0]]->GetBinContent(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu24->SetBinError(1,mu24[binMap[0]]->GetBinError(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu24->SetBinContent(2,mu24[binMap[1]]->GetBinContent(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu24->SetBinError(2,mu24[binMap[1]]->GetBinError(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu24->SetBinContent(3,mu24[binMap[2]]->GetBinContent(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu24->SetBinError(3,mu24[binMap[2]]->GetBinError(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu24->SetBinContent(4,mu24[binMap[3]]->GetBinContent(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu24->SetBinError(4,mu24[binMap[3]]->GetBinError(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu24->SetBinContent(5,mu24[binMap[4]]->GetBinContent(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu24->SetBinError(5,mu24[binMap[4]]->GetBinError(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu24->SetBinContent(6,mu24[binMap[5]]->GetBinContent(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu24->SetBinError(6,mu24[binMap[5]]->GetBinError(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu24->SetBinContent(7,mu24[binMap[6]]->GetBinContent(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu24->SetBinError(7,mu24[binMap[6]]->GetBinError(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu24->SetBinContent(8,mu24[binMap[7]]->GetBinContent(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu24->SetBinError(8,mu24[binMap[7]]->GetBinError(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu24->SetBinContent(9,mu24[binMap[8]]->GetBinContent(1)*scaleFactor_mumu[8]);
  yieldPlot_mumu24->SetBinError(9,mu24[binMap[8]]->GetBinError(1)*scaleFactor_mumu[8]);
  
  TH1D * yieldPlot_mumu = new TH1D("yieldPlot_mumu","",10,0,10);
  yieldPlot_mumu->SetBinContent(1,mu21[binMap[0]]->GetBinContent(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinError(1,mu21[binMap[0]]->GetBinError(1)*scaleFactor_mumu[0]);
  yieldPlot_mumu->SetBinContent(2,mu21[binMap[1]]->GetBinContent(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinError(2,mu21[binMap[1]]->GetBinError(1)*scaleFactor_mumu[1]);
  yieldPlot_mumu->SetBinContent(3,mu21[binMap[2]]->GetBinContent(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinError(3,mu21[binMap[2]]->GetBinError(1)*scaleFactor_mumu[2]);
  yieldPlot_mumu->SetBinContent(4,mu21[binMap[3]]->GetBinContent(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinError(4,mu21[binMap[3]]->GetBinError(1)*scaleFactor_mumu[3]);
  yieldPlot_mumu->SetBinContent(5,mu21[binMap[4]]->GetBinContent(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinError(5,mu21[binMap[4]]->GetBinError(1)*scaleFactor_mumu[4]);
  yieldPlot_mumu->SetBinContent(6,mu21[binMap[5]]->GetBinContent(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinError(6,mu21[binMap[5]]->GetBinError(1)*scaleFactor_mumu[5]);
  yieldPlot_mumu->SetBinContent(7,mu21[binMap[6]]->GetBinContent(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinError(7,mu21[binMap[6]]->GetBinError(1)*scaleFactor_mumu[6]);
  yieldPlot_mumu->SetBinContent(8,mu21[binMap[7]]->GetBinContent(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinError(8,mu21[binMap[7]]->GetBinError(1)*scaleFactor_mumu[7]);
  yieldPlot_mumu->SetBinContent(9,mu21[binMap[8]]->GetBinContent(1)*scaleFactor_mumu[8]);
  yieldPlot_mumu->SetBinError(9,mu21[binMap[8]]->GetBinError(1)*scaleFactor_mumu[8]);
  
  TH1D * yieldPlot_ee = new TH1D("yieldPlot_ee","",10,-0.25,9.75);
  yieldPlot_ee->SetBinContent(1,e[binMap[0]]->GetBinContent(1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinError(1,e[binMap[0]]->GetBinError(1)*scaleFactor_ee[0]);
  yieldPlot_ee->SetBinContent(2,e[binMap[1]]->GetBinContent(1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinError(2,e[binMap[1]]->GetBinError(1)*scaleFactor_ee[1]);
  yieldPlot_ee->SetBinContent(3,e[binMap[2]]->GetBinContent(1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinError(3,e[binMap[2]]->GetBinError(1)*scaleFactor_ee[2]);
  yieldPlot_ee->SetBinContent(4,e[binMap[3]]->GetBinContent(1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinError(4,e[binMap[3]]->GetBinError(1)*scaleFactor_ee[3]);
  yieldPlot_ee->SetBinContent(5,e[binMap[4]]->GetBinContent(1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinError(5,e[binMap[4]]->GetBinError(1)*scaleFactor_ee[4]);
  yieldPlot_ee->SetBinContent(6,e[binMap[5]]->GetBinContent(1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinError(6,e[binMap[5]]->GetBinError(1)*scaleFactor_ee[5]);
  yieldPlot_ee->SetBinContent(7,e[binMap[6]]->GetBinContent(1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinError(7,e[binMap[6]]->GetBinError(1)*scaleFactor_ee[6]);
  yieldPlot_ee->SetBinContent(8,e[binMap[7]]->GetBinContent(1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinError(8,e[binMap[7]]->GetBinError(1)*scaleFactor_ee[7]);
  yieldPlot_ee->SetBinContent(9,e[binMap[8]]->GetBinContent(1)*scaleFactor_ee[8]);
  yieldPlot_ee->SetBinError(9,e[binMap[8]]->GetBinError(1)*scaleFactor_ee[8]);

  TH1D * yieldCombo = new TH1D("yieldCombo","",10,0.25,10.25);
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){//subtract 1 for the pp bin
    yieldCombo->SetBinContent(i, combo[binMap[i-1]]->GetBinContent(1) * scaleFactor_Combo[i-1] ); 
    yieldCombo->SetBinError(i, combo[binMap[i-1]]->GetBinError(1) * scaleFactor_Combo[i-1] );
  }

  TH1D * ATLAS = new TH1D("ATLAS","",1,9,10);
  ATLAS->SetBinContent(1,0.3745*TMath::Power(10,-6));
  ATLAS->SetBinError(1,0.0034*TMath::Power(10,-6));
  ATLAS->SetMarkerStyle(8);
  ATLAS->SetMarkerColor(kGreen+1);
  ATLAS->SetLineColor(kGreen+1);

  yieldPlot_mumu24->SetMarkerStyle(8);
  yieldPlot_mumu24->SetMarkerColor(kBlue+1);
  yieldPlot_mumu24->SetMarkerSize(1.3);
  yieldPlot_mumu24->SetLineColor(kBlue+1);
  yieldPlot_mumu24->SetLineWidth(2);
  
  yieldPlot_mumu->SetMarkerStyle(8);
  yieldPlot_mumu->SetMarkerColor(kViolet+1);
  yieldPlot_mumu->SetMarkerSize(1.3);
  yieldPlot_mumu->SetLineColor(kViolet+1);
  yieldPlot_mumu->SetLineWidth(2);
  yieldPlot_mumu->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z} (mb)");
  yieldPlot_mumu->GetYaxis()->SetRangeUser(0.15*TMath::Power(10,-6),0.45*TMath::Power(10,-6));
  if(doAccept) yieldPlot_mumu->GetYaxis()->SetRangeUser(0,0.45*TMath::Power(10,-6)*1.6);

  yieldPlot_ee->SetMarkerStyle(21);
  yieldPlot_ee->SetMarkerSize(1.5);
  yieldPlot_ee->SetMarkerColor(kRed+1);
  yieldPlot_ee->SetLineColor(kRed+1);
  yieldPlot_ee->SetLineWidth(2);
  
  yieldCombo->SetMarkerStyle(34);
  yieldCombo->SetMarkerSize(1.8);
  yieldCombo->SetMarkerColor(kBlack);
  yieldCombo->SetLineColor(kBlack);
  yieldCombo->SetLineWidth(2);

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  yieldPlot_mumu->SetStats(0);
  
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0)");
  setex2->Draw();

  yieldPlot_mumu->Draw();
  yieldPlot_mumu->Draw("same");
  //yieldPlot_mumu24->Draw("same");
  yieldPlot_ee->Draw("same");
  yieldCombo->Draw("same");
  if(!doAccept) ATLAS->Draw("same"); 
 
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");
  setex1->Draw();  
  
  //HGPythia scaled by 0-90 yield
  TH1D * hgp = new TH1D("HG_PYTHIA","",10,0,10);
  for(int i = 1; i<9; i++){
    float inclusivePbPb = yieldCombo->GetBinContent(9);
    float inclusivePbPb_e = TMath::Sqrt(yieldCombo->GetBinError(9) * yieldCombo->GetBinError(9) + TMath::Power( comboSyst[binMap[9-1]]->GetBinContent(1)* scaleFactor_Combo[9-1] , 2) + TMath::Power( yieldCombo->GetBinContent(9) * TAARelErr[9-1]  ,2));
    hgp->SetBinContent(i,hgPythia->GetBinContent(i) * inclusivePbPb );
    hgp->SetBinError(i,hgPythia->GetBinContent(i) * inclusivePbPb_e );
  }
  hgp->SetFillColor(kGreen-6);
  hgp->SetLineWidth(0);
  hgp->Draw("same E2");
  setex2->Draw();

  TBox * eBox[40];
  TBox * mu21Box[40];
  //TBox * mu24Box[40];
  TBox * netBox[40];
  TBox * ATLASBox;  

  TBox * glauberBox[40][4];
  //TBox * glauberBox2[40][4];

  //*********GLAUBER BRACKETS
  
  float Gwidth = 0.1;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){
      helper.drawBoxAbsolute(yieldCombo, i , glauberBox[i][0], yieldCombo->GetBinContent(i) * TAARelErr[i-1],Gwidth,(Color_t)kGray+1, true, 1001); 
      helper.drawBoxAbsolute(yieldPlot_ee, i , glauberBox[i][1], yieldPlot_ee->GetBinContent(i) * TAARelErr[i-1] ,Gwidth ,(Color_t)kGray+1, true, 1001); 
      helper.drawBoxAbsolute(yieldPlot_mumu, i , glauberBox[i][2], yieldPlot_mumu->GetBinContent(i) * TAARelErr[i-1],Gwidth,(Color_t)kGray+1, true, 1001); 
     // helper.drawBoxAbsolute(yieldPlot_mumu24, i , glauberBox[i][3], yieldPlot_mumu24->GetBinContent(i) * TAARelErr[i-1],Gwidth,(Color_t)kGray+2); 
  }
  /*float shift = 0.005*TMath::Power(10,-6);
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){
      helper.drawBoxAbsolute(yieldCombo, i , glauberBox2[i][0], yieldCombo->GetBinContent(i) * TAARelErr[i-1] - shift,Gwidth,(Color_t)kWhite); 
      helper.drawBoxAbsolute(yieldPlot_ee, i , glauberBox2[i][1], yieldPlot_ee->GetBinContent(i) * TAARelErr[i-1] - shift,Gwidth ,(Color_t)kWhite); 
      helper.drawBoxAbsolute(yieldPlot_mumu, i , glauberBox2[i][2], yieldPlot_mumu->GetBinContent(i) * TAARelErr[i-1]- shift,Gwidth,(Color_t)kWhite); 
    //  helper.drawBoxAbsolute(yieldPlot_mumu24, i , glauberBox2[i][3], yieldPlot_mumu24->GetBinContent(i) * TAARelErr[i-1] - shift,Gwidth,(Color_t)kWhite); 
  }*/
  //**********END GLAUBER BRACKETS
  
  //Syst Uncertainties
  float width = 0.06;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){
      helper.drawBoxAbsolute(yieldCombo, i , netBox[i], comboSyst[binMap[i-1]]->GetBinContent(1)* scaleFactor_Combo[i-1],width,(Color_t)kBlack); 
      helper.drawBoxAbsolute(yieldPlot_ee, i , eBox[i], yieldPlot_ee->GetBinContent(i) * totalError[binMap[i-1]][0]->GetBinContent(1) ,width ,(Color_t)kRed+1); 
      helper.drawBoxAbsolute(yieldPlot_mumu, i , mu21Box[i], yieldPlot_mumu->GetBinContent(i) * totalError[binMap[i-1]][1]->GetBinContent(1),width,(Color_t)kViolet+1); 
    //  helper.drawBoxAbsolute(yieldPlot_mumu24, i , mu24Box[i], yieldPlot_mumu24->GetBinContent(i) * totalError[binMap[i-1]][2]->GetBinContent(1),width,(Color_t)kBlue); 
  }
  

  for(int i = 1; i<10; i++){
    yieldPlot_mumu->GetXaxis()->SetBinLabel(i, labels[i-1]);
    yieldPlot_mumu->GetXaxis()->ChangeLabel(i,45);
  }
  
 
  TH1D * glauberDummy = new TH1D("glauberDummy","",1,0,1);
  glauberDummy->SetLineColor(kGray+1);
  glauberDummy->SetFillColor(kGray+1);
 
  //ATLAS point error    
  if(!doAccept) helper.drawBoxAbsolute(ATLAS, 1 , ATLASBox, 0.00787*TMath::Power(10,-6) ,width,(Color_t)kGreen+1); 
  
  yieldPlot_mumu->Draw("same");
  //yieldPlot_mumu24->Draw("same");
  yieldPlot_ee->Draw("same");
  yieldCombo->Draw("same");
  setex1->Draw();  
 
  yieldPlot_mumu->Print("All");
  yieldPlot_ee->Print("All");
  yieldCombo->Print("All");  

  TLegend * leg = new TLegend(0.21,0.22,0.83,0.50);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"p_{T}^{l} > 20 GeV, |y_{Z}| < 2.1","");
  if(!doAccept){
    //leg->AddEntry(yieldPlot_mumu24,"#mu^{+}#mu^{-} |#eta^{l}| < 2.4","p");
    leg->AddEntry(yieldPlot_mumu,"#mu^{+}#mu^{-} |#eta^{l}| < 2.1","p");
    leg->AddEntry(yieldPlot_ee,"e^{+}e^{-} |#eta^{l}| < 2.1 ","p");
    leg->AddEntry(yieldCombo,"Combined |#eta^{l}| < 2.1","p");
  } else {
    //leg->AddEntry(yieldPlot_mumu24,"#mu^{+}#mu^{-} |y_{Z}| < 2.4","p");
    leg->AddEntry(yieldPlot_mumu,"#mu^{+}#mu^{-}","p");
    leg->AddEntry(yieldPlot_ee,"e^{+}e^{-}","p");
    leg->AddEntry(yieldCombo,"Combined","p");
  }
  leg->AddEntry(glauberDummy,"Glauber Uncertainties","f");
  if(!doAccept) leg->AddEntry(ATLAS,"ATLAS pp |#eta^{l}| < 2.5","p");
  leg->AddEntry(hgp,"Scaled HG-PYTHIA","f");

  leg->Draw("same");
  c1->RedrawAxis();
  CMS_lumi(c1,0,10);
  
  if(!doAccept){
    c1->SaveAs("plots/prettyPlots/yields_Pretty.png");
    c1->SaveAs("plots/prettyPlots/yields_Pretty.pdf");
    c1->SaveAs("plots/prettyPlots/yields_Pretty.C");
  }else{
    c1->SaveAs("plots/prettyPlots/yields_Pretty_withAccept.png");
    c1->SaveAs("plots/prettyPlots/yields_Pretty_withAccept.pdf");
    c1->SaveAs("plots/prettyPlots/yields_Pretty_withAccept.C");
  }

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 7)
  {
    std::cout << "Usage: massPeakPlots <Z2EE file> <Z2mumu21 file> <Z2mumu24 file> < Z2EE syst> <Z2mumu21 syst> <Z2mumu24 syst>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu21 = argv[2];
  std::string Zmumu24 = argv[3];
  std::string ZeeSyst = argv[4];
  std::string Zmumu21Syst = argv[5];
  std::string Zmumu24Syst = argv[6];
   
  plotMassPeaks(Zee, Zmumu21, Zmumu24, ZeeSyst, Zmumu21Syst, Zmumu24Syst, false);
  plotMassPeaks(Zee, Zmumu21, Zmumu24, ZeeSyst, Zmumu21Syst, Zmumu24Syst, true);
  return 0; 
}

#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
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
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

float getSigma(float x){
  if(x>1 || x<0){
    std::cout << "x should be [0,1]" << std::endl;
    return -1;
  }
  for(int i = 0; i<100000; i++){
    float p = 0.5 * (1+TMath::Erf((float)i/10000.0 / TMath::Sqrt(2))); //probability of lying between -infinity and i/10000 sigma;
    float pComplement = 1-p; //probability of lying between i/10000sigma and infinity
    float pTails = 2*pComplement;//probability in the area of both + and - tails
    if( pTails < x ) return (float)i/10000.0;
  }
  return 10;
}

void plotMassPeaks(std::string Zee, std::string Zmumu21, std::string Zmumu24, std::string ZeeSyst, std::string Zmumu21Syst, std::string Zmumu24Syst, bool doAccept){
	  Settings s = Settings();

	  float unitScale = TMath::Power(10,6);
	  //float sigmaNN = 0.487585;
	  //float sigmaNN = 0.43789;


	  CentralityTool c = CentralityTool();
	  const int nBins = c.getNCentBins();

	  CombinePoints cp = CombinePoints();
	  HistNameHelper helper = HistNameHelper();

  TFile * HGPythia = TFile::Open("resources/HGPythia_hfSum_doZ0.root","read");
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

  std::vector< double > eCWeight;
  std::vector< double > muCWeight;

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
      bool isUsed = j<6 || j==15 || j==16;
      if( i==1 && isUsed)  muCWeight.push_back(combined.at(3));
      if( i==1 && isUsed)  eCWeight.push_back(combined.at(4));
    }
  }

  std::cout << eCWeight.size() << " " << muCWeight.size() << std::endl;
  for( unsigned int i = 0; i<eCWeight.size(); i++){
    std::cout << eCWeight.at(i) << " " << muCWeight.at(i) << std::endl;
  }

  const char * labels[10] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%", "50-70%", "70-90%", "0-90%","pp"};
  const char * labelse[10] = {"ee 0-5%","ee 5-10%","ee 10-20%","ee 20-30%","ee 30-40%","ee 40-50%", "ee 50-70%", "ee 70-90%", "ee 0-90%","pp"};
  const char * labelsmu[10] = {"#mu#mu 0-5%","#mu#mu 5-10%","#mu#mu 10-20%","#mu#mu 20-30%","#mu#mu 30-40%","#mu#mu 40-50%", "#mu#mu 50-70%", "#mu#mu 70-90%", "#mu#mu 0-90%","pp"};
  float TAA[9] = {25.70, 20.40, 14.39, 8.798, 5.124, 2.777, 0.9957, 0.1650, 6.274};
  float TAARelErr[9] = {0.018, 0.020, 0.021, 0.025, 0.031, 0.038, 0.050, 0.047, 0.022};
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
 
  yieldPlot_mumu->GetYaxis()->CenterTitle(); 
  yieldPlot_mumu->GetYaxis()->SetTitleSize(0.05); 
  yieldPlot_mumu->GetYaxis()->SetTitleOffset(1.6); 
  yieldPlot_mumu->GetYaxis()->SetLabelSize(0.05); 
  yieldPlot_mumu->GetYaxis()->SetLabelOffset(0.005); 
  yieldPlot_mumu->GetXaxis()->SetTitleSize(0.06); 
  yieldPlot_mumu->GetXaxis()->SetTitleOffset(1.1); 
  yieldPlot_mumu->GetXaxis()->SetLabelSize(0.05); 
  yieldPlot_mumu->GetXaxis()->SetLabelOffset(0.008); 
  yieldPlot_mumu->GetYaxis()->CenterTitle(); 
  yieldPlot_mumu->GetXaxis()->CenterTitle(); 
  yieldPlot_mumu->GetXaxis()->SetTitleOffset(1.1);
  yieldPlot_mumu->GetXaxis()->SetTitle("Centrality"); 
  yieldPlot_mumu->SetMarkerStyle(8);
  yieldPlot_mumu->SetMarkerColor(kBlue);
  yieldPlot_mumu->SetMarkerSize(1.3);
  yieldPlot_mumu->SetLineColor(kBlue);
  yieldPlot_mumu->SetLineWidth(2);
  yieldPlot_mumu->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z} (nb)");

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

  yieldPlot_mumu->Scale(unitScale);
  yieldPlot_ee->Scale(unitScale);
  yieldCombo->Scale(unitScale);
  yieldPlot_mumu->GetYaxis()->SetRangeUser(0.15*TMath::Power(10,-6)*unitScale,0.45*TMath::Power(10,-6)*unitScale*1.6);
  if(doAccept) yieldPlot_mumu->GetYaxis()->SetRangeUser(0,0.72);
  
  yieldPlot_mumu->Draw();

  //uncertainty on theory line
  float theoryLine = yieldCombo->GetBinContent(9);
  float theoryUncert = TMath::Sqrt(yieldCombo->GetBinError(9)*yieldCombo->GetBinError(9) + theoryLine*theoryLine* TAARelErr[9-1]*TAARelErr[9-1] + TMath::Power(comboSyst[binMap[9-1]]->GetBinContent(1)* scaleFactor_Combo[9-1]*unitScale,2) ) ; 
  TBox * sigmaBox = new TBox(yieldPlot_mumu->GetXaxis()->GetBinLowEdge(1),theoryLine+theoryUncert,yieldPlot_mumu->GetXaxis()->GetBinUpEdge(yieldPlot_mumu->GetSize()-2),theoryLine-theoryUncert);
  sigmaBox->SetLineWidth(0);
  sigmaBox->SetFillColorAlpha(kMagenta+2,0.2);
  TLine * line1;

  //theory line
  line1 = new TLine(yieldPlot_mumu->GetXaxis()->GetBinLowEdge(1),theoryLine,yieldPlot_mumu->GetXaxis()->GetBinUpEdge(yieldPlot_mumu->GetSize()-2),theoryLine);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);
  line1->SetLineColor(kMagenta+2);
  sigmaBox->Draw("same");
  line1->Draw("same");
  //TLatex latex;
  //latex.SetTextSize(0.037);
  //latex.DrawLatex(9,0.51,"#sigma^{Z}_{NN}");

  yieldPlot_mumu->Draw("same");
  //yieldPlot_mumu24->Draw("same");
  yieldPlot_ee->Draw("same");
  yieldCombo->Draw("same");
  yieldCombo->Print("All");
  if(!doAccept) ATLAS->Draw("same"); 
 
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");
  setex1->Draw();  
  
  //HGPythia scaled by 0-90 yield
  TH1D * hgp = new TH1D("HG_PYTHIA","",10,0,10);
  TH1D * hgpDiv = new TH1D("HG_PYTHIA_Div","",10,0,10);
  for(int i = 1; i<10; i++){
    //float inclusivePbPb = yieldCombo->GetBinContent(9);
    //float inclusivePbPb_e = TMath::Sqrt(yieldCombo->GetBinError(9) * yieldCombo->GetBinError(9) + TMath::Power( comboSyst[binMap[9-1]]->GetBinContent(1)* scaleFactor_Combo[9-1] , 2) + TMath::Power( yieldCombo->GetBinContent(9) * TAARelErr[9-1]  ,2));
    //float inclusivePbPb_e = TMath::Power(10,-7) * 0.025*unitScale;
    //hgp->SetBinContent(i,hgPythia->GetBinContent(i) * inclusivePbPb );
    //hgp->SetBinError(i,hgPythia->GetBinContent(i) * inclusivePbPb_e );
    hgp->SetBinContent(i,hgPythia->GetBinContent(i) * theoryLine );
    hgpDiv->SetBinContent(i,hgPythia->GetBinContent(i) * theoryLine );
    hgp->SetBinError(i,hgPythia->GetBinContent(i) * (theoryUncert) );
    hgpDiv->SetBinError(i,0 );
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
  std::cout << "Drawing syst boxes" << std::endl;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-2; i++){
      helper.drawBoxAbsolute(yieldCombo, i , netBox[i], comboSyst[binMap[i-1]]->GetBinContent(1)* scaleFactor_Combo[i-1]*unitScale,width,(Color_t)kBlack);
      std::cout << i << " " <<  comboSyst[binMap[i-1]]->GetBinContent(1)* scaleFactor_Combo[i-1]*unitScale/yieldCombo->GetBinContent(i) << std::endl;
      helper.drawBoxAbsolute(yieldPlot_ee, i , eBox[i], yieldPlot_ee->GetBinContent(i) * totalError[binMap[i-1]][0]->GetBinContent(1) ,width ,(Color_t)kRed+1); 
      helper.drawBoxAbsolute(yieldPlot_mumu, i , mu21Box[i], yieldPlot_mumu->GetBinContent(i) * totalError[binMap[i-1]][1]->GetBinContent(1),width,(Color_t)kBlue); 
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
  //setex1->Draw();  
 
  yieldPlot_mumu->Print("All");
  yieldPlot_ee->Print("All");
  yieldCombo->Print("All");  

  TLegend * leg = new TLegend(0.21,0.22,0.83,0.55);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"|y_{Z}| < 2.1","");
  if(!doAccept){
    //leg->AddEntry(yieldPlot_mumu24,"#mu^{+}#mu^{-} |#eta^{l}| < 2.4","p");
    leg->AddEntry(yieldPlot_mumu,"#mu^{+}#mu^{-} |#eta^{l}| < 2.1","p");
    leg->AddEntry(yieldPlot_ee,"e^{+}e^{-} |#eta^{l}| < 2.1 ","p");
    leg->AddEntry(yieldCombo,"Combined |#eta^{l}| < 2.1","p");
  } else {
    //leg->AddEntry(yieldPlot_mumu24,"#mu^{+}#mu^{-} |y_{Z}| < 2.4","p");
    leg->AddEntry(yieldPlot_mumu,"Z/#gamma*#rightarrow #mu^{+}#mu^{-}","p");
    leg->AddEntry(yieldPlot_ee,"Z/gamma*#rightarrow e^{+}e^{-}","p");
    leg->AddEntry(yieldCombo,"Combined","p");
  }
  
  std::cout << "yield Compatibility" << std::endl;
  /*float maxDeviationPT = -99;
  for(int i = 1; i<yieldPlot_mumu->GetSize()-1; i++){
    float e = yieldPlot_ee->GetBinContent(i);
    float mu = yieldPlot_mumu->GetBinContent(i);
    float difference = TMath::Abs(e-mu);
 
    float eStat = yieldPlot_ee->GetBinError(i);
    float muStat = yieldPlot_mumu->GetBinError(i);
    float eSyst = yieldPlot_ee->GetBinContent(i) * totalError[binMap[i-1]][0]->GetBinContent(1);
    float muSyst = yieldPlot_mumu->GetBinContent(i) * totalError[binMap[i-1]][1]->GetBinContent(1);
    float err = TMath::Sqrt( eStat * eStat + muStat * muStat + eSyst * eSyst + muSyst * muSyst );
    std::cout << i << " " << difference/err << " sigma" << std::endl;
      
    if(difference/err > maxDeviationPT) maxDeviationPT = difference/err;
  }*/
  //std::cout << "Max deviation: " << maxDeviationPT << std::endl;

  leg->AddEntry(glauberDummy,"Glauber Uncertainties","f");
  //leg->AddEntry(line1,"#sigma^{Z}_{NN} MG5_aMC@NLO + CT14 + EPPS16","lf");
  leg->AddEntry(line1,"0--90% data","lf");
  if(!doAccept) leg->AddEntry(ATLAS,"ATLAS pp |#eta^{l}| < 2.5","p");
  leg->AddEntry(hgp,"T_{AA}-scaled HG-PYTHIA","f");

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

  //calculating probabilities
  

  float fullyCorr = TMath::Sqrt(theoryUncert*theoryUncert/theoryLine/theoryLine + 0.006*0.006);//0.8% for model, 0.6 for acceptance, 1.1 for Nmb

  //add in the Nmb uncertainty to TAARelErr
  for(int i = 0; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
    if(i<7) TAARelErr[i] = TMath::Sqrt(TAARelErr[i]*TAARelErr[i]+ 0.00611*0.00611);
    if(i==7) TAARelErr[i] = TMath::Sqrt(TAARelErr[i]*TAARelErr[i]+ 0.05*0.05);
    if(i==8) TAARelErr[i] = TMath::Sqrt(TAARelErr[i]*TAARelErr[i]+ 0.0126*0.0126);
  }

  float chi2Flat = 0;
  float chi2HGPythia = 0;
  int ndof = 0;
  float chi2Flat_last3 = 0;
  float chi2HGPythia_last3 = 0;
  int ndof_last3 = 0;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldCombo->GetBinContent(i);
      float statUncert = yieldCombo->GetBinError(i);
      float systUncert =  comboSyst[binMap[i-1]]->GetBinContent(1)* scaleFactor_Combo[i-1]*unitScale;
      float TAAUncert = yieldCombo->GetBinContent(i) * TAARelErr[i-1];
      float totalUncert = TMath::Sqrt( statUncert * statUncert + systUncert * systUncert + TAAUncert * TAAUncert + fullyCorr*fullyCorr*TMath::Power(yieldCombo->GetBinContent(i) - 2*TAARelErr[i-1]*TAARelErr[9-1],2) );//last factor accounts for TAA correlations
      std::cout << i << " " << totalUncert << std::endl;
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);

      if(i==8){
        std::cout << "70-90 percent deviation: " << point-theoryFlat << " " << totalUncert << " "  << (point-theoryFlat)/totalUncert  << std::endl;
        std::cout << "70-90 percent deviation: " << TMath::Power((point-theoryFlat)/totalUncert,2) << " " << TMath::Prob(TMath::Power((point-theoryFlat)/totalUncert,2),1) << " " << getSigma(TMath::Prob(TMath::Power((point-theoryFlat)/totalUncert,2),1)  )  << std::endl;
        std::cout << "70-90 percent deviation HGPYTHIA: " << point-theoryHGPythia << " " << totalUncert << " "  << (point-theoryHGPythia)/totalUncert  << std::endl;
      }

      chi2Flat += TMath::Power( (point-theoryFlat)/totalUncert , 2);
      chi2HGPythia += TMath::Power( (point-theoryHGPythia)/totalUncert, 2);
      ndof++;
      if(i>5){
        chi2Flat_last3 += TMath::Power( (point-theoryFlat)/totalUncert , 2);
        chi2HGPythia_last3 += TMath::Power( (point-theoryHGPythia)/totalUncert, 2);
        ndof_last3++;
      }
  }
  std::cout << " " << std::endl;
  std::cout << "All uncorrelated, combined data" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat << "  ndof - " << ndof << "  chi2/ndof - " << chi2Flat/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat,ndof)<< "  sigma - " << getSigma(TMath::Prob(chi2Flat,ndof)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia << "  ndof - " << ndof << "  chi2/ndof - " << chi2HGPythia/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia,ndof)) << std::endl;
  std::cout << " " << std::endl;
  std::cout << "All uncorrelated, combined data (last3)" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat_last3 << "  ndof - " << ndof_last3 << "  chi2/ndof - " << chi2Flat_last3/(float)ndof_last3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat_last3,ndof_last3)<< "  sigma - " << getSigma(TMath::Prob(chi2Flat_last3,ndof_last3)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia_last3 << "  ndof - " << ndof_last3 << "  chi2/ndof - " << chi2HGPythia_last3/(float)ndof_last3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia_last3,ndof_last3) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia_last3,ndof_last3)) << std::endl;
  

  //with correlations....
  chi2Flat = 0;
  chi2HGPythia = 0;
  TMatrixD covarStat = TMatrixD(ndof,ndof);
  TMatrixD covarSyst = TMatrixD(ndof,ndof);
  TMatrixD covarTAA = TMatrixD(ndof,ndof);
  TMatrixD obsFlat = TMatrixD(ndof,1);
  TMatrixD obsHGPythia = TMatrixD(ndof,1);
  TMatrixD obsFlatT = TMatrixD(1,ndof);
  TMatrixD obsHGPythiaT = TMatrixD(1,ndof);
  ndof = 0;
  //std::vector< float > obsFlat;
  //std::vector< float > obsHGPythia;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldCombo->GetBinContent(i);
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);
      obsFlat[i-1][0] = (point-theoryFlat);
      obsHGPythia[i-1][0] = (point-theoryHGPythia);
      obsFlatT[0][i-1] = (point-theoryFlat);
      obsHGPythiaT[0][i-1] = (point-theoryHGPythia);

      float statUncert = yieldCombo->GetBinError(i);
      float systUncert =  comboSyst[binMap[i-1]]->GetBinContent(1)* scaleFactor_Combo[i-1]*unitScale;
      float TAAUncert = yieldCombo->GetBinContent(i) * TMath::Sqrt(TAARelErr[i-1]*TAARelErr[i-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      
      covarStat[i-1][i-1] = statUncert*statUncert;
      for(int j = 1; j<yieldCombo->GetXaxis()->GetNbins()+2-3; j++){
        if(i==j) covarSyst[i-1][j-1] = systUncert * comboSyst[binMap[j-1]]->GetBinContent(1)* scaleFactor_Combo[j-1]*unitScale;
        covarTAA[i-1][j-1] = TAAUncert * yieldCombo->GetBinContent(j) * TMath::Sqrt( TAARelErr[j-1] * TAARelErr[j-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      }
      ndof++;
  }
  TMatrixD inverseCovar =  covarStat + covarSyst + covarTAA;
  for(int i = 0; i<8; i++){
    for(int j = 0; j<8; j++){
      std::cout << inverseCovar[i][j] << " ";
    }
    std::cout << std::endl;
  }
  inverseCovar.Invert();
  chi2Flat = (float)((obsFlatT * inverseCovar * obsFlat)[0][0]);
  chi2HGPythia = (float)((obsHGPythiaT * inverseCovar * obsHGPythia)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst uncorrelated/TAA etc. Correlated, combined data" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat << "  ndof - " << ndof << "  chi2/ndof - " << chi2Flat/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2Flat,ndof)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia << "  ndof - " << ndof << "  chi2/ndof - " << chi2HGPythia/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia,ndof)) << std::endl;
   
   //last 3 points...
   TMatrixD obsFlat2_last3_2 = TMatrixD(3,1);
   TMatrixD obsHGPythia2_last3_2 = TMatrixD(3,1);
   TMatrixD obsFlatT2_last3_2 = TMatrixD(1,3);
   TMatrixD obsHGPythiaT2_last3_2 = TMatrixD(1,3);
   TMatrixD inverseCovar2_last3_2 = TMatrixD(3,3);
   for(int i = 0; i<3; i++){
     obsFlat2_last3_2[i][0] = obsFlat[i+5][0];
     obsFlatT2_last3_2[0][i] = obsFlatT[0][i+5];
     obsHGPythia2_last3_2[i][0] = obsHGPythia[i+5][0];
     obsHGPythiaT2_last3_2[0][i] = obsHGPythiaT[0][i+5];

     for(int j = 0; j<3; j++){
       inverseCovar2_last3_2[i][j] = covarStat[i+5][j+5] + covarSyst[i+5][j+5] + covarTAA[i+5][j+5]; 
     }
   }
  inverseCovar2_last3_2.Invert();
  float chi2Flat_last3_2 = (float)((obsFlatT2_last3_2 * inverseCovar2_last3_2 * obsFlat2_last3_2)[0][0]);
  float chi2HGPythia_last3_2 = (float)((obsHGPythiaT2_last3_2 * inverseCovar2_last3_2 * obsHGPythia2_last3_2)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst uncorrelated/TAA Correlated, using combined (last 3 points only)" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat_last3_2 << "  ndof - " << 3 << "  chi2/ndof - " << chi2Flat_last3_2/(float)3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat_last3_2,3) << "  sigma - " << getSigma(TMath::Prob(chi2Flat_last3_2,3)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia_last3_2 << "  ndof - " << 3 << "  chi2/ndof - " << chi2HGPythia_last3_2/(float)3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia_last3_2,3) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia_last3_2,3)) << std::endl;
  
  
  //with both channels instead
  chi2Flat = 0;
  chi2HGPythia = 0;
  ndof = 0;
  //ee
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldPlot_ee->GetBinContent(i);
      float statUncert = yieldPlot_ee->GetBinError(i);
      float systUncert =  yieldPlot_ee->GetBinContent(i) * totalError[binMap[i-1]][0]->GetBinContent(1) ;
      float TAAUncert = yieldPlot_ee->GetBinContent(i) * TMath::Sqrt(TAARelErr[i-1] * TAARelErr[i-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      float totalUncert = TMath::Sqrt( statUncert * statUncert + systUncert * systUncert + TAAUncert * TAAUncert );
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);
      std::cout << i << " " << totalUncert << std::endl;
      chi2Flat += TMath::Power( (point-theoryFlat)/totalUncert , 2);
      chi2HGPythia += TMath::Power( (point-theoryHGPythia)/totalUncert, 2);
      ndof++;
  }
  //mumu
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldPlot_mumu->GetBinContent(i);
      float statUncert = yieldPlot_mumu->GetBinError(i);
      float systUncert =  yieldPlot_mumu->GetBinContent(i) * totalError[binMap[i-1]][1]->GetBinContent(1) ;
      float TAAUncert = yieldPlot_mumu->GetBinContent(i) * TMath::Sqrt( TAARelErr[i-1] * TAARelErr[i-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      float totalUncert = TMath::Sqrt( statUncert * statUncert + systUncert * systUncert + TAAUncert * TAAUncert );
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);
      std::cout << i+8 << " " << totalUncert << std::endl;
      chi2Flat += TMath::Power( (point-theoryFlat)/totalUncert , 2);
      chi2HGPythia += TMath::Power( (point-theoryHGPythia)/totalUncert, 2);
      ndof++;
  }
  std::cout << " " << std::endl;
  std::cout << "All uncorrelated, using each separate channel" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat << "  ndof - " << ndof << "  chi2/ndof - " << chi2Flat/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat,ndof)<< "  sigma - " << getSigma(TMath::Prob(chi2Flat,ndof)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia << "  ndof - " << ndof << "  chi2/ndof - " << chi2HGPythia/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia,ndof)) << std::endl;
 

  //best estimation is here 
  chi2Flat = 0;
  chi2HGPythia = 0;
  TMatrixD covarStat2 = TMatrixD(ndof,ndof);
  TMatrixD covarSyst2 = TMatrixD(ndof,ndof);
  TMatrixD covarSyst2_uncorr = TMatrixD(ndof,ndof);
  TMatrixD covarTAA2 = TMatrixD(ndof,ndof);
  TMatrixD obsFlat2 = TMatrixD(ndof,1);
  TMatrixD obsHGPythia2 = TMatrixD(ndof,1);
  TMatrixD obsFlatT2 = TMatrixD(1,ndof);
  TMatrixD obsHGPythiaT2 = TMatrixD(1,ndof);
  int tempNDOF = ndof/2;
  ndof = 0;
  //std::vector< float > obsFlat;
  //std::vector< float > obsHGPythia;
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldPlot_ee->GetBinContent(i);
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);
      obsFlat2[i-1][0] = (point-theoryFlat);
      obsHGPythia2[i-1][0] = (point-theoryHGPythia);
      obsFlatT2[0][i-1] = (point-theoryFlat);
      obsHGPythiaT2[0][i-1] = (point-theoryHGPythia);

      float statUncert = yieldPlot_ee->GetBinError(i);
      float systUncert =  yieldPlot_ee->GetBinContent(i) * totalError[binMap[i-1]][0]->GetBinContent(1);
      float TAAUncert = yieldPlot_ee->GetBinContent(i) * TMath::Sqrt(TAARelErr[i-1] * TAARelErr[i-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      
      //statistical diagonal
      covarStat2[i-1][i-1] = statUncert*statUncert;
      //em background
      float EMuncert = emError[binMap[i-1]][0]->GetBinContent(1) * yieldPlot_ee->GetBinContent(i); 
      covarStat2[i-1][i-1] = covarStat2[i-1][i-1] + EMuncert*EMuncert; 

      //hfUncerti
      float HFuncert_i = hfError[binMap[i-1]][0]->GetBinContent(1);
      
      float NmbUncert = 0;
      if(i==8) NmbUncert = 0.05;
      else if(i==9) NmbUncert = 0.01261;
      else NmbUncert = 0.00611;
 
      HFuncert_i = TMath::Sqrt(HFuncert_i * HFuncert_i - NmbUncert*NmbUncert);//remove the Nmb value
      if(i<4) HFuncert_i *= -1;

      //electron correlation matrix
      TFile * eleCorrf = TFile::Open("resources/Z2ee_EfficiencyMC_0_CorrelationMatrix.root","read");
      TH2D * eleCorr = (TH2D*) eleCorrf->Get("correlation_TnP");
      for(int j = 1; j<yieldCombo->GetXaxis()->GetNbins()+2-3; j++){
        //uncorrelated assumption
        if(i==j) covarSyst2_uncorr[i-1][j-1] = systUncert * yieldPlot_ee->GetBinContent(j) * totalError[binMap[j-1]][0]->GetBinContent(1);

        //correlated assumption, upper left corner syst
        //covarSyst2[i-1][j-1] = systUncert * yieldPlot_ee->GetBinContent(j) * totalError[binMap[j-1]][0]->GetBinContent(1);

        //leptons
        covarSyst2[i-1][j-1] = covarSyst2[i-1][j-1] + efficiencyError[binMap[i-1]][0]->GetBinContent(1) * yieldPlot_ee->GetBinContent(i) * efficiencyError[binMap[j-1]][0]->GetBinContent(1) * yieldPlot_ee->GetBinContent(j) * eleCorr->GetBinContent(i,j);
        
        //charge swapping
        covarSyst2[i-1][j-1] = covarSyst2[i-1][j-1] + 0.005 * yieldPlot_ee->GetBinContent(i) * 0.005 * yieldPlot_ee->GetBinContent(j);

        //hfUncertj
        float HFuncert_j = hfError[binMap[j-1]][0]->GetBinContent(1);
        if(j==8) NmbUncert = 0.05;
        else if(j==9) NmbUncert = 0.01261;
        else NmbUncert = 0.00611;
        HFuncert_j = TMath::Sqrt(HFuncert_j * HFuncert_j - NmbUncert*NmbUncert);//remove the Nmb value
        if(j<4) HFuncert_j *= -1;
        covarSyst2[i-1][j-1] = covarSyst2[i-1][j-1] + HFuncert_i * HFuncert_j * yieldPlot_ee->GetBinContent(i) * yieldPlot_ee->GetBinContent(j);
       
        //extra for the e-mu cross correlation for TAA
        float HFuncert_muj = hfError[binMap[j-1]][1]->GetBinContent(1);
        HFuncert_muj = TMath::Sqrt(HFuncert_muj * HFuncert_muj - NmbUncert*NmbUncert);//remove the Nmb value
        if(j<4) HFuncert_muj *= -1;
        covarSyst2[i-1][j-1+tempNDOF] = covarSyst2[i-1][j-1+tempNDOF] + HFuncert_i * yieldPlot_ee->GetBinContent(i) * HFuncert_muj * yieldPlot_mumu->GetBinContent(j);
        covarSyst2[j-1+tempNDOF][i-1] = covarSyst2[i-1][j-1+tempNDOF];

        //TAA
        covarTAA2[i-1][j-1] = TAAUncert * yieldPlot_ee->GetBinContent(j) * TMath::Sqrt( TAARelErr[j-1] * TAARelErr[j-1] + fullyCorr*fullyCorr);//0.8% is model uncert
        //cross terms
        covarTAA2[i-1][j-1+tempNDOF] = TAAUncert * yieldPlot_mumu->GetBinContent(j) * TMath::Sqrt( TAARelErr[j-1] * TAARelErr[j-1] + fullyCorr*fullyCorr);//0.8% is model uncert
        covarTAA2[j-1+tempNDOF][i-1] = covarTAA2[i-1][j-1+tempNDOF];//other side by symmetry
      }
      ndof++;
  }
  //muon corrleation matrix
  TFile * muCorrf = TFile::Open("resources/Z2mumu_Efficiencies_CorrelationMatrix.root","read");
  TH2D * muCorr = (TH2D*) muCorrf->Get("correlation_TnP");
  for(int i = 1; i<yieldCombo->GetXaxis()->GetNbins()+2-3; i++){
      float point = yieldPlot_mumu->GetBinContent(i);
      float theoryFlat = theoryLine;
      float theoryHGPythia = hgp->GetBinContent(i);
      obsFlat2[i-1+tempNDOF][0] = (point-theoryFlat);
      obsHGPythia2[i-1+tempNDOF][0] = (point-theoryHGPythia);
      obsFlatT2[0][i-1+tempNDOF] = (point-theoryFlat);
      obsHGPythiaT2[0][i-1+tempNDOF] = (point-theoryHGPythia);

      float statUncert = yieldPlot_mumu->GetBinError(i);
      float systUncert =  yieldPlot_mumu->GetBinContent(i) * totalError[binMap[i-1]][1]->GetBinContent(1);
      float TAAUncert = yieldPlot_mumu->GetBinContent(i) * TMath::Sqrt( TAARelErr[i-1] * TAARelErr[i-1] + fullyCorr*fullyCorr - 2*TAARelErr[9-1]*TAARelErr[i-1]);//removing correlated TAA component with minus sign
     
      //statistical diagnoal 
      covarStat2[i-1+tempNDOF][i-1+tempNDOF] = statUncert*statUncert;
      //em background
      float EMuncert = emError[binMap[i-1]][1]->GetBinContent(1) * yieldPlot_mumu->GetBinContent(i); 
      covarStat2[i-1+tempNDOF][i-1+tempNDOF] = covarStat2[i-1+tempNDOF][i-1+tempNDOF] + EMuncert*EMuncert; 
      
      //hfUncerti
      float HFuncert_i = hfError[binMap[i-1]][1]->GetBinContent(1);

      float NmbUncert = 0;
      if(i==8) NmbUncert = 0.05;
      else if(i==9) NmbUncert = 0.01261;
      else NmbUncert = 0.00611;

      HFuncert_i = TMath::Sqrt(HFuncert_i * HFuncert_i - NmbUncert * NmbUncert);//remove the Nmb value
      if(i<4) HFuncert_i *= -1;
      for(int j = 1; j<yieldCombo->GetXaxis()->GetNbins()+2-3; j++){
        //uncorrelated assumption
        if(i==j) covarSyst2_uncorr[i-1+tempNDOF][j-1+tempNDOF] = systUncert * yieldPlot_mumu->GetBinContent(j) * totalError[binMap[j-1]][1]->GetBinContent(1);
   
        //correlated assumption, bottom right corner syst
        //covarSyst2[i-1+tempNDOF][j-1+tempNDOF] = systUncert * yieldPlot_mumu->GetBinContent(j) * totalError[binMap[j-1]][1]->GetBinContent(1);
        //leptons
        covarSyst2[i-1+tempNDOF][j-1+tempNDOF] = covarSyst2[i-1+tempNDOF][j-1+tempNDOF] + efficiencyError[binMap[i-1]][1]->GetBinContent(1) * yieldPlot_mumu->GetBinContent(i) * efficiencyError[binMap[j-1]][1]->GetBinContent(1) * yieldPlot_mumu->GetBinContent(j) * muCorr->GetBinContent(i,j);
        
        //hfUncertj
        float HFuncert_j = hfError[binMap[j-1]][1]->GetBinContent(1);
        if(j==8) NmbUncert = 0.05;
        else if(j==9) NmbUncert = 0.01261;
        else NmbUncert = 0.00611;
        HFuncert_j = TMath::Sqrt(HFuncert_j * HFuncert_j - NmbUncert * NmbUncert);//remove the Nmb value
        if(j<4) HFuncert_j *= -1;
        covarSyst2[i-1+tempNDOF][j-1+tempNDOF] = covarSyst2[i-1+tempNDOF][j-1+tempNDOF] + HFuncert_i * HFuncert_j * yieldPlot_mumu->GetBinContent(i) * yieldPlot_mumu->GetBinContent(j);
   
        //TAA
        covarTAA2[i-1+tempNDOF][j-1+tempNDOF] = TAAUncert * yieldPlot_mumu->GetBinContent(j) * TMath::Sqrt( TAARelErr[j-1] * TAARelErr[j-1] - 2*TAARelErr[9-1]*TAARelErr[j-1] + fullyCorr*fullyCorr);//removing correlated fraction from 0-90% TAA with minus sign
        //extra for the e-mu cross correlation for TAA (already done in the previous iteration from symmetry)
        //covarTAA2[i-1+tempNDOF][j-1] = TAAUncert * yieldPlot_ee->GetBinContent(j) * TMath::Sqrt(TAARelErr[j-1] * TAARelErr[j-1] + fullyCorr*fullyCorr);//0.8% is model uncert
      }
      ndof++;
  }

  //correlation matrix
  TH2D * correlationMatrix = new TH2D("correlationMatrix","",16,0,16,16,0,16);
  for(int i = 1; i<17; i++){
    for(int j = 1; j<17; j++){
      correlationMatrix->SetBinContent(i,j, covarStat2[i-1][j-1] + covarSyst2[i-1][j-1] + covarTAA2[i-1][j-1]);
    }
  }
  float diagonals[16] = {0};
  for( int i = 0; i<16; i++){
     diagonals[i] = TMath::Sqrt(correlationMatrix->GetBinContent(i+1,i+1));
     std::cout << i << " " << diagonals[i] << std::endl;
  }
  for(int i = 1; i<17; i++){
    for(int j = 1; j<17; j++){
      correlationMatrix->SetBinContent(i,j, (covarStat2[i-1][j-1] + covarSyst2[i-1][j-1] + covarTAA2[i-1][j-1]) / diagonals[i-1] / diagonals[j-1]);
    }
  }
  TCanvas * corr = new TCanvas("corr","corr",800,800);
  correlationMatrix->SetStats(0);
  correlationMatrix->GetXaxis()->SetTitle("Bin");
  correlationMatrix->GetYaxis()->SetTitle("Bin");
  for(int i = 1; i<9; i++){
    correlationMatrix->GetXaxis()->SetBinLabel(i, labelse[i-1]);
    correlationMatrix->GetYaxis()->SetBinLabel(i, labelse[i-1]);
    correlationMatrix->GetXaxis()->SetBinLabel(i+8, labelsmu[i-1]);
    correlationMatrix->GetYaxis()->SetBinLabel(i+8, labelsmu[i-1]);
  }
  correlationMatrix->Draw("colz");
  corr->SaveAs("plots/correlationMatrix.png");
  corr->SaveAs("plots/correlationMatrix.pdf");

  TMatrixD inverseCovar2 =  covarStat2 + covarSyst2 + covarTAA2;
  TMatrixD inverseCovar2Copy = inverseCovar2;
  inverseCovar2.Invert();
  chi2Flat = (float)((obsFlatT2 * inverseCovar2 * obsFlat2)[0][0]);
  chi2HGPythia = (float)((obsHGPythiaT2 * inverseCovar2 * obsHGPythia2)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst/TAA Correlated, using each separate channel" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat << "  ndof - " << ndof << "  chi2/ndof - " << chi2Flat/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2Flat,ndof)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia << "  ndof - " << ndof << "  chi2/ndof - " << chi2HGPythia/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia,ndof)) << std::endl;
 
   //last 3 points...
   TMatrixD obsFlat2_last3 = TMatrixD(6,1);
   TMatrixD obsHGPythia2_last3 = TMatrixD(6,1);
   TMatrixD obsFlatT2_last3 = TMatrixD(1,6);
   TMatrixD obsHGPythiaT2_last3 = TMatrixD(1,6);
   TMatrixD inverseCovar2_last3 = TMatrixD(6,6);
   for(int i = 0; i<3; i++){
     obsFlat2_last3[i][0] = obsFlat2[i+5][0];
     obsFlat2_last3[i+3][0] = obsFlat2[i+13][0];
     obsFlatT2_last3[0][i] = obsFlatT2[0][i+5];
     obsFlatT2_last3[0][i+3] = obsFlatT2[0][i+13];
     obsHGPythia2_last3[i][0] = obsHGPythia2[i+5][0];
     obsHGPythia2_last3[i+3][0] = obsHGPythia2[i+13][0];
     obsHGPythiaT2_last3[0][i] = obsHGPythiaT2[0][i+5];
     obsHGPythiaT2_last3[0][i+3] = obsHGPythiaT2[0][i+13];

     for(int j = 0; j<3; j++){
       inverseCovar2_last3[i][j] = inverseCovar2Copy[i+5][j+5]; 
       inverseCovar2_last3[i+3][j+3] = inverseCovar2Copy[i+13][j+13]; 
       inverseCovar2_last3[i][j+3] = inverseCovar2Copy[i+5][j+13]; 
       inverseCovar2_last3[i+3][j] = inverseCovar2Copy[i+13][j+5]; 
     }
   }
  inverseCovar2_last3.Invert();
  chi2Flat_last3 = (float)((obsFlatT2_last3 * inverseCovar2_last3 * obsFlat2_last3)[0][0]);
  chi2HGPythia_last3 = (float)((obsHGPythiaT2_last3 * inverseCovar2_last3 * obsHGPythia2_last3)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst/TAA Correlated, using each separate channel (last 3 points only)" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat_last3 << "  ndof - " << 6 << "  chi2/ndof - " << chi2Flat_last3/(float)6 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat_last3,6) << "  sigma - " << getSigma(TMath::Prob(chi2Flat_last3,6)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia_last3 << "  ndof - " << 6 << "  chi2/ndof - " << chi2HGPythia_last3/(float)6 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia_last3,6) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia_last3,6)) << std::endl;
   
   //combined points from full corr matrix...
   std::cout << "starting combo with full corr" << std::endl;
   std::cout << eCWeight.size() << " " << muCWeight.size() << std::endl;
   TMatrixD ComboCovar2 = TMatrixD(8,8);
   for(int i = 0; i<8; i++){
     for(int j = 0; j<8; j++){
       ComboCovar2[i][j] = inverseCovar2Copy[i][j] * eCWeight.at(i) * eCWeight.at(j) + inverseCovar2Copy[i+8][j] * muCWeight.at(i) * eCWeight.at(j) + inverseCovar2Copy[i][j+8] * eCWeight.at(i) * muCWeight.at(j) + inverseCovar2Copy[i+8][j+8] * muCWeight.at(i) * muCWeight.at(j);
       if(i==j) std::cout << i << " " << TMath::Sqrt(ComboCovar2[i][j]) << std::endl;
     }
   }
  for(int i = 0; i<8; i++){
    for(int j = 0; j<8; j++){
      std::cout << ComboCovar2[i][j] << " ";
    }
    std::cout << std::endl;
  }
   TMatrixD inverseComboCovar2 = ComboCovar2;
   inverseComboCovar2.Invert();
   std::cout << "inversion done" << std::endl;  
 
   TMatrixD ComboFlat2 = TMatrixD(8,1);
   TMatrixD ComboHGPythia2 = TMatrixD(8,1);
   TMatrixD ComboFlatT2 = TMatrixD(1,8);
   TMatrixD ComboHGPythiaT2 = TMatrixD(1,8);
  for(int i = 1; i<9; i++){
    float point = yieldCombo->GetBinContent(i);
    float theoryFlat = theoryLine;
    float theoryHGPythia = hgp->GetBinContent(i);
    ComboFlat2[i-1][0] = point-theoryFlat;
    ComboHGPythia2[i-1][0] =  point-theoryHGPythia;
    ComboFlatT2[0][i-1] = point-theoryFlat;
    ComboHGPythiaT2[0][i-1] =  point-theoryHGPythia;
  }
  float Combochi2Flat = (float)((ComboFlatT2 * inverseComboCovar2 * ComboFlat2)[0][0]);
  float Combochi2HGPythia = (float)((ComboHGPythiaT2 * inverseComboCovar2 * ComboHGPythia2)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst/TAA Correlated, using combined channel but full covar matrix" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << Combochi2Flat << "  ndof - " << 8 << "  chi2/ndof - " << Combochi2Flat/(float)8 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(Combochi2Flat,8) << "  sigma - " << getSigma(TMath::Prob(Combochi2Flat,8)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << Combochi2HGPythia << "  ndof - " << 8 << "  chi2/ndof - " << Combochi2HGPythia/(float)8 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(Combochi2HGPythia,8) << "  sigma - " << getSigma(TMath::Prob(Combochi2HGPythia,8)) << std::endl;

  //last 3
  TMatrixD Combo3 = TMatrixD(3,3);
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      Combo3[i][j] = ComboCovar2[i+5][j+5]; 
    }
  }
  TMatrixD inverseCombo3 = Combo3;
  inverseCombo3.Invert();
  
   TMatrixD Combo3Flat2 = TMatrixD(3,1);
   TMatrixD Combo3HGPythia2 = TMatrixD(3,1);
   TMatrixD Combo3FlatT2 = TMatrixD(1,3);
   TMatrixD Combo3HGPythiaT2 = TMatrixD(1,3);
  for(int i = 0; i<3; i++){
    float point = yieldCombo->GetBinContent(i+6);
    std::cout << point << std::endl;
    float theoryFlat = theoryLine;
    float theoryHGPythia = hgp->GetBinContent(i+6);
    Combo3Flat2[i][0] = point-theoryFlat;
    Combo3HGPythia2[i][0] =  point-theoryHGPythia;
    Combo3FlatT2[0][i] = point-theoryFlat;
    Combo3HGPythiaT2[0][i] =  point-theoryHGPythia;
  }
  float Combo3chi2Flat = (float)((Combo3FlatT2 * inverseCombo3 * Combo3Flat2)[0][0]);
  float Combo3chi2HGPythia = (float)((Combo3HGPythiaT2 * inverseCombo3 * Combo3HGPythia2)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "Syst/TAA Correlated, using combined channel but full covar matrix (Last 3 points)" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << Combo3chi2Flat << "  ndof - " << 3 << "  chi2/ndof - " << Combo3chi2Flat/(float)3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(Combo3chi2Flat,3) << "  sigma - " << getSigma(TMath::Prob(Combo3chi2Flat,3)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << Combo3chi2HGPythia << "  ndof - " << 3 << "  chi2/ndof - " << Combo3chi2HGPythia/(float)3 << "  Prob (chi2/ndof > observed) - " << TMath::Prob(Combo3chi2HGPythia,3) << "  sigma - " << getSigma(TMath::Prob(Combo3chi2HGPythia,3)) << std::endl;


  //correlation matrix 8x8
  TH2D * correlationMatrix8x8 = new TH2D("correlationMatrix8x8","",8,0,8,8,0,8);
  for(int i = 1; i<9; i++){
    for(int j = 1; j<9; j++){
      correlationMatrix8x8->SetBinContent(i,j, ComboCovar2[i-1][j-1] );
    }
  }
  for( int i = 0; i<8; i++){
     diagonals[i] = TMath::Sqrt(correlationMatrix8x8->GetBinContent(i+1,i+1));
     std::cout << i << " " << diagonals[i] << std::endl;
  }
  for(int i = 1; i<9; i++){
    for(int j = 1; j<9; j++){
      correlationMatrix8x8->SetBinContent(i,j, (ComboCovar2[i-1][j-1]) / diagonals[i-1] / diagonals[j-1]);
    }
  }
  TCanvas * corr8x8 = new TCanvas("corr8x8","corr8x8",800,800);
  correlationMatrix8x8->SetStats(0);
  correlationMatrix8x8->GetXaxis()->SetTitle("Bin");
  correlationMatrix8x8->GetYaxis()->SetTitle("Bin");
  for(int i = 1; i<9; i++){
    correlationMatrix8x8->GetXaxis()->SetBinLabel(i, labels[i-1]);
    correlationMatrix8x8->GetYaxis()->SetBinLabel(i, labels[i-1]);
  }
  correlationMatrix8x8->Draw("colz");
  corr8x8->SaveAs("plots/correlationMatrix8x8.png");
  corr8x8->SaveAs("plots/correlationMatrix8x8.pdf");


   //end of best estimation

  
  TMatrixD inverseCovar2_uncorrSyst =  covarStat2 + covarSyst2_uncorr + covarTAA2;
  inverseCovar2_uncorrSyst.Invert();
  chi2Flat = (float)((obsFlatT2 * inverseCovar2_uncorrSyst * obsFlat2)[0][0]);
  chi2HGPythia = (float)((obsHGPythiaT2 * inverseCovar2_uncorrSyst * obsHGPythia2)[0][0]);
  std::cout << " " << std::endl;
  std::cout << "TAA Correlated, using each separate channel" << std::endl;
  std::cout << "Flat Hypothesis: chi2 - " << chi2Flat << "  ndof - " << ndof << "  chi2/ndof - " << chi2Flat/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2Flat,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2Flat,ndof)) << std::endl;
  std::cout << "HGPythia Hypothesis: chi2 - " << chi2HGPythia << "  ndof - " << ndof << "  chi2/ndof - " << chi2HGPythia/(float)ndof << "  Prob (chi2/ndof > observed) - " << TMath::Prob(chi2HGPythia,ndof) << "  sigma - " << getSigma(TMath::Prob(chi2HGPythia,ndof)) << std::endl;


  //revised plot fo PRL
  
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  c2->SetBorderSize(0);
  TPad * p1 = new TPad("p1","p1",0,0.23,1,1,0);
  TPad * p2 = new TPad("p2","p2",0,0,1,0.23,0);
  c2->SetLineWidth(0);
  p1->SetBottomMargin(0);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.001);
  p1->SetBorderSize(0);
  p1->Draw();
  p2->SetTopMargin(0);
  p2->SetLeftMargin(0.15);
  p2->SetRightMargin(0.001);
  p2->SetBottomMargin(0.5);
  p2->SetBorderSize(0);
  p2->Draw();
  p1->cd();

  const int nGPts = 9;
  float xOff = 0.5;
  float x[nGPts] = {2.5,7.5,15,25,35,45,60,80,100};
  float xErr[nGPts] = {2.5-xOff, 2.5-xOff, 5.-xOff, 5.-xOff, 5.-xOff, 5.-xOff, 10.-xOff, 10.-xOff, 0};
  float xErr_TAA[nGPts] = {2,2,2,2,2,2,2,2,2};
  float xErr_Syst[nGPts] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
  float xErr_hgp[nGPts] = {2.5, 2.5, 5., 5., 5., 5., 10., 10., 7.5};
  float y[nGPts] = {0};
  float y_hgp[nGPts] = {0};
  float y_hgpDiv[nGPts] = {0};
  float yErr[nGPts] = {0};
  float yErr_hgp[nGPts] = {0};
  float yErr_hgpDiv[nGPts] = {0};
  float yErr_hgpDivTAA[nGPts] = {0};
  float yErr_hgpDivSyst[nGPts] = {0};
  float yErr_Syst[nGPts] = {0};
  float yErr_TAA[nGPts] = {0};
  for(int i = 0; i<nGPts; i++){
    y[i] = yieldCombo->GetBinContent(i+1);
    yErr[i] = yieldCombo->GetBinError(i+1);
    y_hgp[i] = hgp->GetBinContent(i+1);
    yErr_hgp[i] = hgp->GetBinError(i+1);
    yErr_TAA[i] = yieldCombo->GetBinContent(i+1) * TAARelErr[i];
    yErr_Syst[i] = comboSyst[binMap[i]]->GetBinContent(1)* scaleFactor_Combo[i]*unitScale ;
  }

  TGraphAsymmErrors * cg = new TGraphAsymmErrors(nGPts, x, y, xErr, xErr, yErr, yErr);  
  TGraphAsymmErrors * cg_TAA = new TGraphAsymmErrors(nGPts,x,y,xErr_TAA,xErr_TAA,yErr_TAA,yErr_TAA);
  TGraphAsymmErrors * cg_Syst = new TGraphAsymmErrors(nGPts,x,y,xErr_Syst,xErr_Syst,yErr_Syst,yErr_Syst);

  TGraphAsymmErrors * hgp_G = new TGraphAsymmErrors(nGPts,x,y_hgp,xErr_hgp,xErr_hgp,yErr_hgp,yErr_hgp);
  cg->SetMarkerStyle(8);
  cg->SetMarkerSize(1.5);
  cg->SetMarkerColor(kBlack);
  cg->SetLineColor(kBlack);
  cg->SetLineWidth(1);  
  cg->GetYaxis()->CenterTitle();
  cg->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{T_{AA}} N_{Z} (nb)");
  cg->GetYaxis()->SetTitleSize(0.07);
  cg->GetYaxis()->SetTitleOffset(0.9);
  cg->GetYaxis()->SetLabelSize(0.055);
  cg->GetYaxis()->SetRangeUser(0.201,0.6);
  cg->GetXaxis()->SetRangeUser(-1,110);
  cg->SetTitle("");
  cg->GetXaxis()->SetTickLength(0);
  cg->Draw("p a");  
  
  sigmaBox->SetX2(111);
  sigmaBox->SetX1(-1);
  sigmaBox->Draw("same");

  hgp_G->SetFillColor(kGreen-6);
  cg_TAA->SetFillColor(kGray+1);
  cg_Syst->SetLineColor(kBlack);
  cg_Syst->SetFillStyle(0);
  gStyle->SetHatchesLineWidth(2);
  hgp_G->SetFillStyle(3344);
  hgp_G->Draw("2 same");
  cg_TAA->Draw("2 same");
  cg_Syst->Draw("5 same");

  line1->SetX2(110);
  line1->SetX1(-1);
  line1->Draw("same");
  TLine * line4;
  line4 = new TLine(90,0.201,90,0.6);
  line4->SetLineWidth(1);
  line4->SetLineStyle(1);
  line4->Draw("same");  
  
  TH1D * dummyLine1 = new TH1D("dummyLine1","",1,0,1);
  dummyLine1->SetFillColorAlpha(kMagenta+2,0.2);
  dummyLine1->SetLineColor(kMagenta+2);
  dummyLine1->SetLineStyle(2);
  dummyLine1->SetLineWidth(2);

  cg->Draw("p same");  

  TLegend * leg3 = new TLegend(0.08,0.7,0.57,0.85);
  leg3->AddEntry((TObject*)0,"60 < m_{ll} < 120 GeV","");
  leg3->AddEntry((TObject*)0,"|y_{Z}| < 2.1","");
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.075);
  leg3->Draw("same");
  
  TLegend * leg2 = new TLegend(0.16,0.04,0.74,0.47);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.075);
  cg->SetFillColor(kGray+1);
  hgp_G->SetLineWidth(0);
  leg2->AddEntry(cg,"Z/#gamma* #rightarrow l^{+}l^{-}","lpef"); 
  //leg2->AddEntry(dummyLine1,"#sigma^{Z}_{NN} MG5_aMC@NLO + CT14 + EPPS16","fl");
  leg2->AddEntry(dummyLine1,"0-90% data","fl");
  leg2->AddEntry(hgp_G,"HG-PYTHIA scaled by 0-90%","f");
  leg2->Draw("same");
 
  p2->cd();
  TH1D * hgpRatio = (TH1D*) yieldCombo->Clone("hgpRatio");
  hgpRatio->Divide(hgpDiv);
  for(int i = 0; i<nGPts; i++){
    y_hgpDiv[i] = hgpRatio->GetBinContent(i+1);
    yErr_hgpDiv[i] = hgpRatio->GetBinError(i+1);
    yErr_hgpDivTAA[i] = hgpRatio->GetBinContent(i+1)* TAARelErr[i];
    yErr_hgpDivSyst[i] = hgpRatio->GetBinContent(i+1) * yErr_Syst[i]/y[i] ;
  }

  TGraphAsymmErrors * hgpRatioG = new TGraphAsymmErrors(nGPts, x, y_hgpDiv, xErr, xErr, yErr_hgpDiv, yErr_hgpDiv);
  TGraphAsymmErrors * hgpRatioGTAA = new TGraphAsymmErrors(nGPts, x, y_hgpDiv, xErr_TAA, xErr_TAA, yErr_hgpDivTAA, yErr_hgpDivTAA);
  TGraphAsymmErrors * hgpRatioGSyst = new TGraphAsymmErrors(nGPts, x, y_hgpDiv, xErr_Syst, xErr_Syst, yErr_hgpDivSyst, yErr_hgpDivSyst);
  hgpRatioG->SetMarkerStyle(8);
  hgpRatioG->SetMarkerSize(1.5);
  hgpRatioG->SetMarkerColor(kBlack);
  hgpRatioG->SetLineColor(kBlack);
  hgpRatioG->SetLineWidth(1);  
  hgpRatioG->GetXaxis()->CenterTitle();  
  hgpRatioG->GetYaxis()->CenterTitle();  
  hgpRatioG->GetYaxis()->SetTitle("#frac{Data}{Model}");  
  hgpRatioG->GetYaxis()->SetTitleSize(0.17);  
  hgpRatioG->GetYaxis()->SetTitleOffset(0.35);  
  hgpRatioG->GetXaxis()->SetTitleSize(0.3);  
  hgpRatioG->GetXaxis()->SetTitleOffset(0.8);  
  hgpRatioG->GetYaxis()->SetLabelSize(0.17);  
  hgpRatioG->GetXaxis()->SetLabelSize(0.2);  
  hgpRatioG->GetYaxis()->SetRangeUser(0.75,1.25);  
  hgpRatioG->GetXaxis()->SetRangeUser(-1,110);  
  hgpRatioG->GetXaxis()->SetTitle("Centrality (%)");  
  hgpRatioG->SetTitle("");  
  hgpRatioG->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
  TAxis* a = hgpRatioG->GetXaxis();
  a->ChangeLabel(6,-1,-1,-1,-1,-1,"0-90");
  hgpRatioG->GetXaxis()->SetTickLength(0);
  hgpRatioG->Draw("p a same");
 
  hgpRatioGTAA->SetFillColor(kGray+1);
  hgpRatioGSyst->SetLineColor(kBlack);
  hgpRatioGSyst->SetFillStyle(0);
  hgpRatioGTAA->Draw("2 same");
  hgpRatioGSyst->Draw("5 same");
 
  hgpRatioG->Draw("p same");

  TLine * line3;
  line3 = new TLine(-1,1.0,110,1.0);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);
  line3->Draw("same");
  TLine * line5 = new TLine(90,0.75,90,1.25);
  line5->SetLineWidth(1);
  line5->SetLineStyle(1);
  line5->Draw("same");  

  TGaxis * aClone4 = new TGaxis(0,0.75,90,0.75,0,90,505,"+");
  aClone4->SetTickLength(0.02);
  aClone4->SetLabelSize(0);
  TGaxis * aClone3 = new TGaxis(0,1.25,90,1.25,0,90,505,"-");
  aClone3->SetTickLength(0.02);
  aClone3->SetLabelSize(0);

  aClone3->Draw();
  aClone4->Draw();
  
  p1->cd();
  TGaxis * aClone = new TGaxis(0,0.2,90,0.2,0,90,505,"+");
  aClone->SetTickLength(0.02);
  aClone->SetLabelSize(0);
  TGaxis * aClone2 = new TGaxis(0,0.6,90,0.6,0,90,505,"-");
  aClone2->SetTickLength(0.02);
  aClone2->SetLabelSize(0);

  aClone->Draw();
  aClone2->Draw();

  c2->cd();
  p1->cd();
  CMS_lumi(p1,0,0,1.5,true,true,false);

  c2->SaveAs("plots/prettyPlots/yields_Pretty_withAccept_PRL.png");
  c2->SaveAs("plots/prettyPlots/yields_Pretty_withAccept_PRL.pdf");
  c2->SaveAs("plots/prettyPlots/yields_Pretty_withAccept_PRL.C");
  
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

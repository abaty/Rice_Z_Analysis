#include "include/HistNameHelper.h"
#include "include/ZEfficiency.h"
#include "include/centralityBin.h"
#include "include/ptReweightSpectrum.h"
#include "include/MCWeightHelper.h"
#include "include/electronEnergyScale.h"
#include "include/forceConsistency.h"
#include "include/electronSelector.h"
#include "include/electronTriggerMatching.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/Timer.h"
#include "include/MCReweight.h"
#include "include/ElectronTnP.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

int pdfProcess(int id1,int id2,int nPartons);

void doZ2EE(std::vector< std::string > files, int jobNumber, bool isTest){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  TH1::SetDefaultSumw2();
  ZEfficiency zEff = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", 0);
  ElectronEnergyScale energyScale = ElectronEnergyScale("MC");
  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronTnP eTnP = ElectronTnP();  

  HistNameHelper histHelper = HistNameHelper();

  MCReweight vzRW = MCReweight("resources/vzReweight.root","resources/centralityFlatteningWeight.root");
  MCWeightHelper weightHelper = MCWeightHelper();
  
  PtReweightSpectrum spectrumRW = PtReweightSpectrum("resources/ptSpectrumReweighting.root");
  
  Settings s = Settings();

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();
  
  TRandom3 * r = new TRandom3();
  TH1D * accept21_pt_pass[120];
  TH1D * accept21_y_pass[120];
  TH1D * accept21_yields_pass[120];
  TH1D * accept21_pt_net[120];
  TH1D * accept21_y_net[120];
  TH1D * accept21_yields_net[120];
  
  TH1D * accept21_pt_noPtWeight_pass;
  TH1D * accept21_y_noPtWeight_pass;
  TH1D * accept21_yields_noPtWeight_pass;
  TH1D * accept21_pt_noPtWeight_net;
  TH1D * accept21_y_noPtWeight_net;
  TH1D * accept21_yields_noPtWeight_net;
  
  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakOS_noSF[nBins]; 
  TH1D * massPeakSS[nBins]; 
  
  TH1D * recoEff_pt_pass[nBins];
  TH1D * recoEff_pt_net[nBins];
  TH1D * recoEff_pt[nBins];  
  TEfficiency * eff_pt[nBins];
  
  TH1D * recoEff_pt_noPtWeight_pass[nBins];
  TH1D * recoEff_pt_noPtWeight_net[nBins];
  TH1D * recoEff_pt_noPtWeight[nBins];  
  TEfficiency * eff_pt_noPtWeight[nBins];
  
  TH1D * recoEff_y_pass[nBins];
  TH1D * recoEff_y_net[nBins];
  TH1D * recoEff_y_net_smeared[nBins];
  TH1D * recoEff_y[nBins];  
  TEfficiency * eff_y[nBins];

  //EE EB BB heck
  TH1D * recoEff_yEE_pass[nBins];
  TH1D * recoEff_yEE_net[nBins];
  TH1D * recoEff_yEE[nBins];  
  TH1D * recoEff_yEB_pass[nBins];
  TH1D * recoEff_yEB_net[nBins];
  TH1D * recoEff_yEB[nBins];  
  TH1D * recoEff_yBB_pass[nBins];
  TH1D * recoEff_yBB_net[nBins];
  TH1D * recoEff_yBB[nBins];  

  
  TH1D * recoEff_phi_pass[nBins];
  TH1D * recoEff_phi_net[nBins];
  TH1D * recoEff_phi[nBins];  
  TEfficiency * eff_phi[nBins];
  
  TH1D * recoEff_cent_pass[nBins];
  TH1D * recoEff_cent_net[nBins];
  TH1D * recoEff_cent[nBins];  
  TEfficiency * eff_cent[nBins];

  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoEff[nBins];
  TEfficiency * eff[nBins];
  
  TH2D * recoEff_U_pass[nBins];
  TH2D * recoEff_U[nBins];
  TEfficiency * eff_U[nBins];
  
  TH2D * recoEff_D_pass[nBins];
  TH2D * recoEff_D[nBins];
  TEfficiency * eff_D[nBins];
  
  TH2D * recoEff_photonU_pass[nBins];
  TH2D * recoEff_photonU[nBins];
  TEfficiency * eff_photonU[nBins];
  
  TH2D * recoEff_photonD_pass[nBins];
  TH2D * recoEff_photonD[nBins];
  TEfficiency * eff_photonD[nBins];
  
  TH2D * recoEff_noSF_pass[nBins];
  TH2D * recoEff_noSF_net[nBins];
  TH2D * recoEff_noSF[nBins];
  TEfficiency * eff_noSF[nBins];

  TH1D * recoEff_pt_pass_forReso_Reco[nBins];
  TH1D * recoEff_pt_pass_forReso_RecoSmeared[nBins];
  TH1D * recoEff_pt_pass_forReso_Gen[nBins];
  TH1D * recoEff_ptReso[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_Reco[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_RecoSmeared[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_NominalToSmeared[nBins];
  
  TH1D * yReso[nBins];
  TH2D * yResponse;
  TH1D * yReco;
  TH1D * yGen;
  TH1D * yRatio;

  TH1D * unfolding_genPt;
  TH1D * unfolding_recoPt;
  TH2D * unfolding_response;
  TH2D * unfolding_response_noPtWeight;
  TH1D * unfolding_genPtU;
  TH1D * unfolding_recoPtU;
  TH2D * unfolding_responseU;
  TH1D * unfolding_genPtD;
  TH1D * unfolding_recoPtD;
  TH2D * unfolding_responseD;

  TH2D * covariance = new TH2D("correlation_TnP","correlation_TnP",8,0,8,8,0,8);

  const int pdfXnbins = 40;
  const double pdfXbins[pdfXnbins] = {0.0001,0.00015,0.0002,0.00025,0.0003,0.0004,0.0005,0.0006,0.0008,0.001,0.003,0.005,0.007,0.01,0.013,0.016,0.02,0.028,0.036,0.044,0.052,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.5, 0.6, 0.7, 0.85, 1.0};

  TH2D * pdfXVsPt = new TH2D("pdfXVsPt","pdfXVsPt",s.nZPtBins-1,s.zPtBins,pdfXnbins-1,pdfXbins);
  TH2D * pdfXVsRap = new TH2D("pdfXVsRap","pdfXVsRap",s.nZRapBins,-s.maxZRap,s.maxZRap,pdfXnbins-1,pdfXbins);
  TH1D * pdfIDVsPt[13];
  for(int i = 0; i<13; i++){
    pdfIDVsPt[i] = new TH1D(Form("pdfIDVsPt_%d",i),Form("pdfIDVsPt_%d",i),s.nZPtBins-1,s.zPtBins);
  }
  TH2D * pdfX1X2LowPt = new TH2D("pdfX1X2LowPt",";x_{projectile};x_{target}",pdfXnbins-1,pdfXbins,pdfXnbins-1,pdfXbins);
  TH2D * pdfX1X2HighPt = new TH2D("pdfX1X2HighPt",";x_{projectile};x_{target}",pdfXnbins-1,pdfXbins,pdfXnbins-1,pdfXbins);

  std::vector<double> theoryEventWeights = std::vector<double>(1210, 0);
  TH1D * theoryPtCT14[57];
  for(int i = 0; i<57; i++) theoryPtCT14[i] = new TH1D(Form("theoryCT14_pt_%d",i),Form("theoryCT14_pt_%d",i),s.nZPtBins-1,s.zPtBins);
  TH1D * theoryPtnCTEQ15[33];
  for(int i = 0; i<33; i++) theoryPtnCTEQ15[i] = new TH1D(Form("theorynCTEQ15_pt_%d",i),Form("theorynCTEQ15_pt_%d",i),s.nZPtBins-1,s.zPtBins);
  TH1D * theoryPtEPPS16[97];
  for(int i = 0; i<97; i++) theoryPtEPPS16[i] = new TH1D(Form("theoryEPPS16_pt_%d",i),Form("theoryEPPS16_pt_%d",i),s.nZPtBins-1,s.zPtBins);
  TH1D * theoryPtScaleVariations[9];
  for(int i = 0; i<9; i++) theoryPtScaleVariations[i] = new TH1D(Form("theoryScaleVariations_pt_%d",i),Form("theoryScaleVariations_pt_%d",i),s.nZPtBins-1,s.zPtBins);
  TH1D * theoryRapCT14[57];
  for(int i = 0; i<57; i++) theoryRapCT14[i] = new TH1D(Form("theoryCT14_rap_%d",i),Form("theoryCT14_rap_%d",i),s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D * theoryRapnCTEQ15[33];
  for(int i = 0; i<33; i++) theoryRapnCTEQ15[i] = new TH1D(Form("theorynCTEQ15_rap_%d",i),Form("theorynCTEQ15_rap_%d",i),s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D * theoryRapEPPS16[97];
  for(int i = 0; i<97; i++) theoryRapEPPS16[i] = new TH1D(Form("theoryEPPS16_rap_%d",i),Form("theoryEPPS16_rap_%d",i),s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D * theoryRapScaleVariations[9];
  for(int i = 0; i<9; i++) theoryRapScaleVariations[i] = new TH1D(Form("theoryScaleVariations_rap_%d",i),Form("theoryScaleVariations_rap_%d",i),s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D * theoryNetEPPS16[97];
  for(int i = 0; i<97; i++) theoryNetEPPS16[i] = new TH1D(Form("theoryEPPS16_Net_%d",i),Form("theoryEPPS16_Net_%d",i),1,0,1);
  TH1D * theoryNetScaleVariations[9];
  for(int i = 0; i<9; i++) theoryNetScaleVariations[i] = new TH1D(Form("theoryScaleVariations_Net_%d",i),Form("theoryScaleVariations_Net_%d",i),1,0,1);

  for(int k = 0; k<nBins; k++){
    massPeakOS[k] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_noSF[k] = new TH1D(Form("massPeakOS_noSF_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[k] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_U_pass[k] = new TH2D(Form("recoEff_U_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_D_pass[k] = new TH2D(Form("recoEff_D_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_photonU_pass[k] = new TH2D(Form("recoEff_photonU_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_photonD_pass[k] = new TH2D(Form("recoEff_photonD_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_noSF_pass[k] = new TH2D(Form("recoEff_noSF_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_noSF_net[k] = new TH2D(Form("recoEff_noSF_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_pt_pass[k] = new TH1D(Form("recoEff_pt_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_pt_net[k] = new TH1D(Form("recoEff_pt_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_pt_noPtWeight_pass[k] = new TH1D(Form("recoEff_pt_noPtWeight_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_pt_noPtWeight_net[k] = new TH1D(Form("recoEff_pt_noPtWeight_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_y_pass[k] = new TH1D(Form("recoEff_y_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_y_net[k] = new TH1D(Form("recoEff_y_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_y_net_smeared[k] = new TH1D(Form("recoEff_y_net_smeared_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_phi_pass[k] = new TH1D(Form("recoEff_phi_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_phi_net[k] = new TH1D(Form("recoEff_phi_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_cent_pass[k] = new TH1D(Form("recoEff_cent_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
    recoEff_cent_net[k] = new TH1D(Form("recoEff_cent_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
 
    //EE EB BB check
    recoEff_yEE_pass[k] = new TH1D(Form("recoEff_yEE_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_yEE_net[k] = new TH1D(Form("recoEff_yEE_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_yEB_pass[k] = new TH1D(Form("recoEff_yEB_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_yEB_net[k] = new TH1D(Form("recoEff_yEB_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_yBB_pass[k] = new TH1D(Form("recoEff_yBB_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_yBB_net[k] = new TH1D(Form("recoEff_yBB_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
 
     
    recoEff_pt_pass_forReso_Reco[k] = new TH1D(Form("recoEff_pt_pass_forReso_Reco_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_RecoSmeared[k] = new TH1D(Form("recoEff_pt_pass_forReso_RecoSmeared_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_Gen[k] = new TH1D(Form("recoEff_pt_pass_forReso_Gen_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_ptReso[k] = new TH1D(Form("recoEff_ptReso_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.2,0.2);
    
    yReso[k] = new TH1D(Form("yReso_%d_%d",c.getCentBinLow(k), c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.1,0.1);
  }
  for(int i = 0; i<weightHelper.getSize(); i++){
    accept21_pt_pass[i] = new TH1D(Form("accept21_pt_pass_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept21_pt_net[i] = new TH1D(Form("accept21_pt_net_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept21_y_pass[i] = new TH1D(Form("accept21_y_pass_%d",i),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
    accept21_y_net[i] = new TH1D(Form("accept21_y_net_%d",i),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
    accept21_yields_pass[i] = new TH1D(Form("accept21_yields_pass_%d",i),"",1,0,1);
    accept21_yields_net[i] = new TH1D(Form("accept21_yields_net_%d",i),"",1,0,1);
  }
  accept21_pt_noPtWeight_pass = new TH1D(Form("accept21_pt_noPtWeight_pass_%d",0),"",s.nZPtBins-1,s.zPtBins);
  accept21_pt_noPtWeight_net = new TH1D(Form("accept21_pt_noPtWeight_net_%d",0),"",s.nZPtBins-1,s.zPtBins);
  accept21_y_noPtWeight_pass = new TH1D(Form("accept21_y_noPtWeight_pass_%d",0),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  accept21_y_noPtWeight_net = new TH1D(Form("accept21_y_noPtWeight_net_%d",0),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  accept21_yields_noPtWeight_pass = new TH1D(Form("accept21_yields_noPtWeight_pass_%d",0),"",1,0,1);
  accept21_yields_noPtWeight_net = new TH1D(Form("accept21_yields_noPtWeight_net_%d",0),"",1,0,1);
  yReco = new TH1D("yReco",";y_{reco}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  yGen = new TH1D("yGen",";y_{gen}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  yResponse = new TH2D("yResponse","y_{reco};y_{gen}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  
  unfolding_genPt = new TH1D("unfolding_genPt","",s.nZPtBins-1,s.zPtBins);
  unfolding_recoPt = new TH1D("unfolding_recoPt","",s.nZPtBins-1,s.zPtBins);
  unfolding_response = new TH2D("unfolding_response","",s.nZPtBins-1,s.zPtBins,s.nZPtBins-1,s.zPtBins);
  unfolding_response_noPtWeight = new TH2D("unfolding_response_noPtWeight","",s.nZPtBins-1,s.zPtBins,s.nZPtBins-1,s.zPtBins);
  unfolding_genPtU = new TH1D("unfolding_genPtU","",s.nZPtBins-1,s.zPtBins);
  unfolding_recoPtU = new TH1D("unfolding_recoPtU","",s.nZPtBins-1,s.zPtBins);
  unfolding_responseU = new TH2D("unfolding_responseU","",s.nZPtBins-1,s.zPtBins,s.nZPtBins-1,s.zPtBins);
  unfolding_genPtD = new TH1D("unfolding_genPtD","",s.nZPtBins-1,s.zPtBins);
  unfolding_recoPtD = new TH1D("unfolding_recoPtD","",s.nZPtBins-1,s.zPtBins);
  unfolding_responseD = new TH2D("unfolding_responseD","",s.nZPtBins-1,s.zPtBins,s.nZPtBins-1,s.zPtBins);

  int nEle;
  int hiBin;
  float hiHF;
  float vz;
  int nMEPartons;
  std::pair<float, float> * pdfX = 0;
  std::pair<int, int> * pdfID = 0;
  std::vector< float > * ttbar_w = 0;

  int pprimaryVertexFilter;
  int phfCoincFilter2Th4;
  int pclusterCompatibilityFilter;

  int HLT_DoubleEle10;
  int HLT_SingleEle20;

  std::vector< float > * elePt = 0;
  std::vector< float > * eleEta = 0;
  std::vector< float > * elePhi = 0;
  std::vector< float > * eleSigmaIEtaIEta = 0;
  std::vector< int > * eleCharge = 0;
  std::vector< int > * eleMissHits = 0;
  //std::vector< float > * eledEtaAtVtx = 0;
  std::vector< float > * eledEtaSeedAtVtx = 0;
  std::vector< float > * eledPhiAtVtx = 0;
  std::vector< float > * eleSCEta = 0;
  std::vector< float > * eleSCPhi = 0;
  std::vector< float > * eleHoverEBc = 0;
  //std::vector< float > * eleD0 = 0;
  //std::vector< float > * eleDz = 0;
  std::vector< float > * eleIP3D = 0;
  std::vector< float > * eleEoverPInv = 0;

  int nMC = 0;
  //std::vector< int > * mcStatus = 0;
  std::vector< int > * mcPID = 0;
  std::vector< float > * mcPt = 0;
  std::vector< float > * mcEta = 0;
  //std::vector< float > * mcPhi = 0;
  std::vector< int > * mcMomPID = 0;
  std::vector< int > * mcGMomPID = 0;
  std::vector< float > * mcMomPt = 0;
  std::vector< float > * mcMomEta = 0;
  std::vector< float > * mcMomPhi = 0; 
  std::vector< float > * mcMomMass = 0;

  unsigned int nFiles = files.size();
  if(isTest) nFiles = 662;
  //if(isTest) nFiles = 20;
  for(unsigned int f = 0; f<nFiles; f++){
    timer.StartSplit("Opening Files");

    TFile * in = TFile::Open(files.at(f).c_str(),"read");
    if(f%5 == 0)  std::cout << f << "/" << files.size() << std::endl;

    TTree * hltTree = (TTree*)in->Get("hltanalysis/HltTree");
    hltTree->SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v1",&HLT_DoubleEle10); 
    hltTree->SetBranchAddress("HLT_HIEle20Gsf_v1",&HLT_SingleEle20); 

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    TTree * eTreeMC = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    eTreeMC->SetBranchAddress("nMC",&nMC);
    //eTreeMC->SetBranchAddress("mcStatus",&mcStatus);
    eTreeMC->SetBranchAddress("mcPID",&mcPID);
    eTreeMC->SetBranchAddress("mcMomPID",&mcMomPID);
    eTreeMC->SetBranchAddress("mcGMomPID",&mcGMomPID);
   
    eTreeMC->SetBranchAddress("mcPt",&mcPt);
    eTreeMC->SetBranchAddress("mcEta",&mcEta);
    //eTreeMC->SetBranchAddress("mcPhi",&mcPhi);
    eTreeMC->SetBranchAddress("mcMomPt",&mcMomPt);
    eTreeMC->SetBranchAddress("mcMomEta",&mcMomEta);
    eTreeMC->SetBranchAddress("mcMomPhi",&mcMomPhi);
    eTreeMC->SetBranchAddress("mcMomMass",&mcMomMass);

    eTree->SetBranchAddress("nEle",&nEle);
    eTree->SetBranchAddress("elePt",&elePt);
    eTree->SetBranchAddress("eleEta",&eleEta);
    eTree->SetBranchAddress("elePhi",&elePhi);
    eTree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta);
    eTree->SetBranchAddress("eleMissHits",&eleMissHits);
    eTree->SetBranchAddress("eleCharge",&eleCharge);
    //eTree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
    eTree->SetBranchAddress("eledEtaSeedAtVtx",&eledEtaSeedAtVtx);
    eTree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
    eTree->SetBranchAddress("eleSCEta",&eleSCEta);
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverEBc",&eleHoverEBc);
    //eTree->SetBranchAddress("eleD0",&eleD0);
    //eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleIP3D",&eleIP3D);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    //evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("hiHF",&hiHF);
    evtTree->SetBranchAddress("vz",&vz);
    evtTree->SetBranchAddress("pdfX",&pdfX);
    evtTree->SetBranchAddress("pdfID",&pdfID);
    evtTree->SetBranchAddress("nMEPartons",&nMEPartons);
    evtTree->SetBranchAddress("ttbar_w",&ttbar_w); 
 
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

    //TTree * L1Tree = (TTree*)in->Get("l1object/L1UpgradeFlatTree");
    //L1Tree->SetBranchAddress("nEGs",&(eTrig.L1nEGs));
    //L1Tree->SetBranchAddress("egEta", &(eTrig.L1egEta));
    //L1Tree->SetBranchAddress("egPhi", &(eTrig.L1egPhi));
    //L1Tree->SetBranchAddress("egEt", &(eTrig.L1egEt));

    TTree * HLTObjTree;
    HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIEle20Gsf_v");
    HLTObjTree->SetBranchAddress("eta",&(eTrig.HLTEta));
    HLTObjTree->SetBranchAddress("phi",&(eTrig.HLTPhi));
    HLTObjTree->SetBranchAddress("pt",&(eTrig.HLTPt));

    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      if(i%5000 == 0) std::cout << i << "/" << eTree->GetEntries() << std::endl;
      timer.StartSplit("Checking Evt Selections");
      //event selections
      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;

      skimTree->GetEntry(i);
      if(! (pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter)) continue;
      
      hiBin = cb.getHiBinFromhiHF(hiHF,3);
      double acceptWeight[120];
      for(int w = 0; w<=1209; w++) theoryEventWeights.at(w) += ttbar_w->at(w);
      for(int w = 0; w<weightHelper.getSize(); w++){
        acceptWeight[w] = vzRW.reweightFactor( vz ) * (ttbar_w->at(weightHelper.getIndx(w))/10000.0);
      }
      double eventWeight = vzRW.reweightFactor( vz ) * vzRW.reweightFactorCent(hiBin) * c.findNcoll( hiBin ) * (ttbar_w->at(1080)/10000.0);//1080 is EEPS16NLO+CT14;
      double eventWeight_noPtWeight = eventWeight;     
 
      timer.StartSplit("Loading GEN electron tree");
      eTreeMC->GetEntry(i);

      timer.StartSplit("Gen Loop (PreElectronConfirmation)");
      int nGenElectronsFound = 0;
      
      TLorentzVector mom = TLorentzVector();
      double ptWeight = 1.0;
      double ptWeightU = 1.0;
      double ptWeightD = 1.0;
      bool foundGen = false;
      bool foundZforTheory = false;

      double firstEleEta = -99;
      double firstElePt = -1;
      for(int j = 0; j<nMC; j++){
        //break out if you find a Z->tautau
        if( mcGMomPID->at(j) == 23 && TMath::Abs(mcMomPID->at(j)) == 15) break;

        if(!foundZforTheory && TMath::Abs(mcMomPID->at(j)) == 23){
          TLorentzVector tempMomThry = TLorentzVector();
          tempMomThry.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
          
          if(TMath::Abs( tempMomThry.Rapidity() ) < s.maxZRap ){
            pdfXVsPt->Fill(mcMomPt->at(j) , pdfX->first, ttbar_w->at(1080));
            pdfXVsRap->Fill(tempMomThry.Rapidity() , pdfX->first, ttbar_w->at(1080));
            pdfIDVsPt[pdfProcess(pdfID->first, pdfID->second, nMEPartons)]->Fill(mcMomPt->at(j), ttbar_w->at(1080));
            pdfIDVsPt[12]->Fill(mcMomPt->at(j), ttbar_w->at(1080));
            if( mcMomPt->at(j) < 10 ) pdfX1X2LowPt->Fill(pdfX->first, pdfX->second, ttbar_w->at(1080));
            if( mcMomPt->at(j) >100 ) pdfX1X2HighPt->Fill(pdfX->first, pdfX->second, ttbar_w->at(1080));

            for(int w = 282; w<=338; w++){
              if(TMath::Abs( tempMomThry.Rapidity() ) < s.maxZRapEle ) theoryPtCT14[w-282]->Fill(mcMomPt->at(j),ttbar_w->at(w));
              theoryRapCT14[w-282]->Fill(tempMomThry.Rapidity(),ttbar_w->at(w));
            }
            for(int w = 1177; w<=1209; w++){
              if(TMath::Abs( tempMomThry.Rapidity() ) < s.maxZRapEle ) theoryPtnCTEQ15[w-1177]->Fill(mcMomPt->at(j),ttbar_w->at(w));
              theoryRapnCTEQ15[w-1177]->Fill(tempMomThry.Rapidity(),ttbar_w->at(w));
            }
            for(int w = 1080; w<=1176; w++){
              if(TMath::Abs( tempMomThry.Rapidity() ) < s.maxZRapEle ){
                theoryPtEPPS16[w-1080]->Fill(mcMomPt->at(j),ttbar_w->at(w));
                theoryNetEPPS16[w-1080]->Fill(0.5,ttbar_w->at(w));
              }
              theoryRapEPPS16[w-1080]->Fill(tempMomThry.Rapidity(),ttbar_w->at(w));
            }
            for(int w = 0; w<=8; w++){
              if(TMath::Abs( tempMomThry.Rapidity() ) < s.maxZRapEle ){
                theoryPtScaleVariations[w]->Fill(mcMomPt->at(j),ttbar_w->at(w));
                theoryNetScaleVariations[w]->Fill(0.5,ttbar_w->at(w));
              }
              theoryRapScaleVariations[w]->Fill(tempMomThry.Rapidity(),ttbar_w->at(w));
            }
          }
          foundZforTheory = true;
        }
        
        //only look for electrons coming directly from a Z
        if( TMath::Abs(mcMomPID->at(j)) != 23) continue;
        
        //if it decays to muons, break out
        if( TMath::Abs(mcPID->at(j)) == 13) break;       
   
        //make sure it's an electron
        if( TMath::Abs(mcPID->at(j)) != 11) continue;
     
        //make sure it's in our rapidity range
        TLorentzVector tempMom = TLorentzVector();
        tempMom.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
        if(TMath::Abs( tempMom.Rapidity()) > s.maxZRapEle ) break;

        //there is a good Z here, let's fill the acceptance histogram
        if(nGenElectronsFound==0){
          double ptWeight = spectrumRW.getReweightFactorElectron(mcMomPt->at(j));
          for(int w = 0; w<weightHelper.getSize(); w++){
            accept21_yields_net[w]->Fill( 0.5, acceptWeight[w]*ptWeight); 
            accept21_y_net[w]->Fill( tempMom.Rapidity(), acceptWeight[w]*ptWeight); 
            accept21_pt_net[w]->Fill( tempMom.Pt(), acceptWeight[w]*ptWeight); 
          } 
          accept21_yields_noPtWeight_net->Fill( 0.5, acceptWeight[0]); 
          accept21_y_noPtWeight_net->Fill( tempMom.Rapidity(), acceptWeight[0]); 
          accept21_pt_noPtWeight_net->Fill( tempMom.Pt(), acceptWeight[0]); 
        }     
        nGenElectronsFound++;
        if(nGenElectronsFound==1){
          firstEleEta = mcEta->at(j);
          firstElePt = mcPt->at(j);
        }

        //if they are not in our acceptance, break out
        if( mcPt->at(j) < s.minElectronPt ) break; 
        if( TMath::Abs(mcEta->at(j)) > s.maxZRapEle ) break; 
      
        timer.StartSplit("Gen Loop (PostElectronConfirmation)");

        //we found both daughters and they are in our acceptance, lets fill our histogram
        if( nGenElectronsFound == 2 ){
          //Fill denominator
          foundGen = true; 
          mom.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
                
          //adjust the event weight based on the gen pT so that the pT spectrum is reweighted to data
          ptWeight = spectrumRW.getReweightFactorElectron(mcMomPt->at(j));
          ptWeightU = spectrumRW.getReweightFactorElectron(mcMomPt->at(j), 1);
          ptWeightD = spectrumRW.getReweightFactorElectron(mcMomPt->at(j), -1);
          eventWeight *= ptWeight;
   
          for(int w = 0; w<weightHelper.getSize(); w++){
            accept21_yields_pass[w]->Fill(0.5, acceptWeight[w] * ptWeight);
            accept21_y_pass[w]->Fill(tempMom.Rapidity(), acceptWeight[w] * ptWeight);
            accept21_pt_pass[w]->Fill(tempMom.Pt(), acceptWeight[w] * ptWeight);
          }
          accept21_yields_noPtWeight_pass->Fill(0.5, acceptWeight[0]);
          accept21_y_noPtWeight_pass->Fill(tempMom.Rapidity(), acceptWeight[0]);
          accept21_pt_noPtWeight_pass->Fill(tempMom.Pt(), acceptWeight[0]);
         
          eTnP.fillCovariance(firstElePt, firstEleEta, mcPt->at(j), mcEta->at(j) ,covariance);
          //if(i%1000==0) covariance->Print("All"); 
          
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              recoEff_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_noSF_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_pt_net[k]->Fill( mom.Pt(), eventWeight);
              recoEff_pt_noPtWeight_net[k]->Fill( mom.Pt(), eventWeight/ptWeight);
              recoEff_y_net[k]->Fill( mom.Rapidity(), eventWeight);
              recoEff_y_net_smeared[k]->Fill( mom.Rapidity() * r->Gaus(1,0.01) , eventWeight);
              //EE EB BB check
              if( TMath::Abs(mcEta->at(j)) < 1.5 && TMath::Abs(firstEleEta) < 1.5 ) recoEff_yBB_net[k]->Fill( mom.Rapidity(), eventWeight);
              else if( TMath::Abs(mcEta->at(j)) > 1.5 && TMath::Abs(firstEleEta) > 1.5 )  recoEff_yEE_net[k]->Fill( mom.Rapidity(), eventWeight);
              else recoEff_yEB_net[k]->Fill( mom.Rapidity(), eventWeight);

              recoEff_cent_net[k]->Fill( hiBin/2.0, eventWeight);
              recoEff_phi_net[k]->Fill( mom.Phi(), eventWeight);

              //for the last few pt bins, halve the y binning so we have better stats
              if(mcMomPt->at(j) > s.zPtBins4Eff[s.nZPtBins4Eff - s.nPtBinsToRebinRapEff]){
                int bin = recoEff_net[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                if( bin%2 ==1){
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_noSF_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                } else {
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_noSF_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
                }
              }
            }
          }
          break; 
        }//end of if statement
      }//end of gen loop      
      if( !foundGen ) continue;
      
      timer.StartSplit("Loading RECO electron tree");
      eTree->GetEntry(i);

      timer.StartSplit("Checking Number of electrons");
      if(nEle<2) continue;
      
      timer.StartSplit("Checking HLT Selections");
      //check for the trigger we want (double ele 10 or single ele 20)
      hltTree->GetEntry(i);
      if( !HLT_SingleEle20) continue;

      timer.StartSplit("Electron Cuts");
      //make a list of electrons passing our cuts
      std::vector< int > goodElectrons;
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        //correct Pt
        elePt->at(j) = energyScale.correctPt( elePt->at(j), eleSCEta->at(j), hiBin);
        
        if(elePt->at(j)< s.minElectronPt) continue;
        if(TMath::Abs(eleSCEta->at(j)) > s.maxZRapEle) continue;
        //veto on transition region
        if(TMath::Abs(eleSCEta->at(j)) > 1.442 && TMath::Abs(eleSCEta->at(j)) < 1.556 ) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        float dEta = TMath::Abs( eledEtaSeedAtVtx->at(j) );
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverEBc->at(j), eleEoverPInv->at(j), eleIP3D->at(j) )) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      bool moreThan2 = false;
      if(goodElectrons.size()>2) moreThan2 = true;

      timer.StartSplit("Loading HLT/L1 Object stuff");
      //get trigger matching stuff
      //L1Tree->GetEntry(i);
      HLTObjTree->GetEntry(i);
    
      //make Z candidates 
      timer.StartSplit("Z candidates");
      TLorentzVector * elec1 = new TLorentzVector();
      TLorentzVector * elec2 = new TLorentzVector();
      for(unsigned int j = 0; j<goodElectrons.size(); j++){
        elec1->SetPtEtaPhiM(elePt->at(goodElectrons.at(j)), eleEta->at(goodElectrons.at(j)), elePhi->at(goodElectrons.at(j)), 0.000511);
        for(unsigned int j2 = j+1; j2<goodElectrons.size(); j2++){
          elec2->SetPtEtaPhiM(elePt->at(goodElectrons.at(j2)), eleEta->at(goodElectrons.at(j2)), elePhi->at(goodElectrons.at(j2)), 0.000511);
          TLorentzVector Zcand = *elec1+*elec2;
          if(Zcand.M() < s.zMassRange[0] || Zcand.M() > s.zMassRange[1]) continue;      

          //acoplanarity cut
          bool passesAco[3] = {1 , 1, 1};
          float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( elePhi->at(goodElectrons.at(j2)) - elePhi->at(goodElectrons.at(j)) )))/TMath::Pi(); 
          if( Zcand.Pt() < s.minPtCutForPhotons_ELE && acoplanarity < s.acoCutForPhotons_ELE ) passesAco[0] = false;
          if( Zcand.Pt() < s.minPtCutForPhotonsU_ELE && acoplanarity < s.acoCutForPhotonsU_ELE ) passesAco[1] = false;
          if( Zcand.Pt() < s.minPtCutForPhotonsD_ELE && acoplanarity < s.acoCutForPhotonsD_ELE ) passesAco[2] = false;
          
          //L1 trigger matching (1 L1 EG > 15 GeV) (remove this requirement on Aug 6th after talking with EGamma POG
          //bool isFirstElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 15.0);
          //bool isSecondElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 15.0);
          //if(! (isFirstElectronL1Matched || isSecondElectronL1Matched)) continue;

          //HLT trigger matching (1 HLT match > 20 GeV)
          bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
          bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
          if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
  

          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " " << Zcand.Eta() << " " << Zcand.Phi() << " " << mom.Pt() << " " << mom.Eta() << " " << mom.Phi() <<  " isOS? " << (int)isOppositeSign << std::endl;
          

           //Looks like this is a good candidate match! Let's get the scale factor
          float scaleFactor = eTnP.getZSF(hiBin, elePt->at(goodElectrons.at(j)), eleSCEta->at(goodElectrons.at(j)), elePt->at(goodElectrons.at(j2)), eleSCEta->at(goodElectrons.at(j2)), 0) ;
          float scaleFactorU = eTnP.getZSF(hiBin, elePt->at(goodElectrons.at(j)), eleSCEta->at(goodElectrons.at(j)), elePt->at(goodElectrons.at(j2)), eleSCEta->at(goodElectrons.at(j2)), 1) ;
          float scaleFactorD = eTnP.getZSF(hiBin, elePt->at(goodElectrons.at(j)), eleSCEta->at(goodElectrons.at(j)), elePt->at(goodElectrons.at(j2)), eleSCEta->at(goodElectrons.at(j2)), -1) ;
          if(moreThan2 && TMath::ACos(TMath::Cos(Zcand.Phi() - mom.Phi())) > 0.1 ) continue;
          if(!isOppositeSign){
            for(int k = 0; k<nBins; k++){ 
              massPeakSS[k]->Fill( Zcand.M(), eventWeight*scaleFactor );
            }
            continue;
          }

          //Fill numerator (and apply the scale factor here!)
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
              if( mom.Pt() < s.zPtBins4Eff[ s.nZPtBins4Eff-1 ] && TMath::Abs( mom.Rapidity() ) < s.maxZRap ){
                if(passesAco[0]){

                  //90%
                  if(k==25){
                    yReco->Fill( Zcand.Rapidity(), eventWeight);
                    yGen->Fill( mom.Rapidity(), eventWeight);
                    yResponse->Fill(Zcand.Rapidity(), mom.Rapidity(), eventWeight);
                  }


                  massPeakOS[k]->Fill(Zcand.M(), eventWeight*scaleFactor);
                  massPeakOS_noSF[k]->Fill(Zcand.M(), eventWeight);
                  recoEff_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                  recoEff_noSF_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                  recoEff_U_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactorU );
                  recoEff_D_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactorD );
                  recoEff_pt_pass[k]->Fill( mom.Pt(), eventWeight * scaleFactor);
                  recoEff_pt_noPtWeight_pass[k]->Fill(  mom.Pt(), eventWeight/ptWeight * scaleFactor );
                  recoEff_y_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                  
                  //EE EB BB check
                  if(TMath::Abs( eleEta->at(goodElectrons.at(j)) ) < 1.5 && TMath::Abs( eleEta->at(goodElectrons.at(j2)) ) < 1.5 ) recoEff_yBB_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                  else if(TMath::Abs( eleEta->at(goodElectrons.at(j)) ) > 1.5 && TMath::Abs( eleEta->at(goodElectrons.at(j2)) ) > 1.5) recoEff_yEE_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                  else recoEff_yEB_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                  
                  recoEff_cent_pass[k]->Fill( hiBin/2.0, eventWeight * scaleFactor);
                  recoEff_phi_pass[k]->Fill( mom.Phi(), eventWeight * scaleFactor);
             
                  recoEff_pt_pass_forReso_Reco[k]->Fill(Zcand.Pt(), eventWeight);         
                  recoEff_pt_pass_forReso_RecoSmeared[k]->Fill(Zcand.Pt() * r->Gaus(1,0.05), eventWeight);         
                  recoEff_pt_pass_forReso_Gen[k]->Fill(mom.Pt(), eventWeight);         
                  recoEff_ptReso[k]->Fill( (Zcand.Pt() - mom.Pt()) / mom.Pt(), eventWeight);         
                
                  yReso[k]->Fill( Zcand.Rapidity() - mom.Rapidity() , eventWeight);            
                
                  if(k==25){//only fill on 0-90% selection
                    float efficiency = zEff.getEfficiency(mom.Rapidity(), mom.Pt(), hiBin);
                    unfolding_recoPt->Fill(Zcand.Pt(),eventWeight*scaleFactor/efficiency);
                    unfolding_genPt->Fill(mom.Pt(),eventWeight*scaleFactor/efficiency);
                    unfolding_response->Fill(Zcand.Pt(), mom.Pt() ,eventWeight*scaleFactor/efficiency);
                    unfolding_response_noPtWeight->Fill(Zcand.Pt(), mom.Pt() ,eventWeight_noPtWeight*scaleFactor/efficiency);
 
                    if(mom.Pt() > Zcand.Pt() + 80){
                      std::cout << "Response anomaly detected!" << std::endl;
                      std::cout << mom.Pt() << " " << mom.Eta() << " " << mom.Phi() << " " << Zcand.Pt() << " " << Zcand.Eta() << " " << Zcand.Phi() << std::endl;
                      std::cout << f << " " <<  i << " " << j << " " << j2 << std::endl;
                    }
                  
                    unfolding_recoPtU->Fill(Zcand.Pt(),eventWeight*scaleFactor*ptWeightU/ptWeight/efficiency);
                    unfolding_genPtU->Fill(mom.Pt(),eventWeight*scaleFactor*ptWeightU/ptWeight/efficiency);
                    unfolding_responseU->Fill(Zcand.Pt(), mom.Pt() ,eventWeight*scaleFactor*ptWeightU/ptWeight/efficiency);
                    unfolding_recoPtD->Fill(Zcand.Pt(),eventWeight*scaleFactor*ptWeightD/ptWeight/efficiency);
                    unfolding_genPtD->Fill(mom.Pt(),eventWeight*scaleFactor*ptWeightD/ptWeight/efficiency);
                    unfolding_responseD->Fill(Zcand.Pt(), mom.Pt() ,eventWeight*scaleFactor*ptWeightD/ptWeight/efficiency);
                  }
          
                  //for the last few pt bins, halve the y binning so we have better stats
                  if(mom.Pt()> s.zPtBins4Eff[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                    int bin = recoEff_pass[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                    if( bin%2 ==1){
                      recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactor );
                      recoEff_U_pass[k]->Fill( recoEff_U_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactorU );
                      recoEff_D_pass[k]->Fill( recoEff_U_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactorD );
                      recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight );
                    } else { 
                      recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactor );
                      recoEff_U_pass[k]->Fill( recoEff_U_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactorU );
                      recoEff_D_pass[k]->Fill( recoEff_D_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactorD );
                      recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight );
                    }
                  }
                }//aco if statement
                if(passesAco[1]){
                  recoEff_photonU_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                  if(mom.Pt()> s.zPtBins4Eff[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                    int bin = recoEff_pass[k]->GetXaxis()->FindBin(mom.Rapidity()); 
                    if( bin%2 ==1){
                      recoEff_photonU_pass[k]->Fill( recoEff_photonU_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactor );
                    }else{
                      recoEff_photonU_pass[k]->Fill( recoEff_photonU_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactor );
                    }
                  }
                }
                if(passesAco[2]){
                  recoEff_photonD_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                  if(mom.Pt()> s.zPtBins4Eff[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                    int bin = recoEff_pass[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                    if( bin%2 ==1){
                      recoEff_photonD_pass[k]->Fill( recoEff_photonD_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactor );
                    }else{
                      recoEff_photonD_pass[k]->Fill( recoEff_photonD_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactor );
                    }
                  }
                }
              }
            }
          }
        }
      }
      delete elec1;
      delete elec2;
    }

    delete eTree;
    in->Close();
  }

  timer.StartSplit("End of analysis");
  std::vector< bool > isConsistent;
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_pass[i], recoEff_net[i]);
    recoEff[i] = (TH2D*)recoEff_pass[i]->Clone(Form("recoEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff[i]->Divide(recoEff_net[i]);
    recoEff[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pass[i]), *(recoEff_net[i]),"w") ){
      eff[i] = new TEfficiency(*(recoEff_pass[i]), *(recoEff_net[i]));
      eff[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff[i]->SetName(Form("eff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff[i]->SetDirectory(0);
      isConsistent.push_back(true);
    }
    else{
      isConsistent.push_back(false);
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }
    recoEff_net[i]->SetDirectory(0);
    recoEff_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_U_pass[i], recoEff_net[i]);
    recoEff_U[i] = (TH2D*)recoEff_U_pass[i]->Clone(Form("recoEff_U_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_U[i]->Divide(recoEff_net[i]);
    recoEff_U[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_U_pass[i]), *(recoEff_net[i]),"w") ){
      eff_U[i] = new TEfficiency(*(recoEff_U_pass[i]), *(recoEff_net[i]));
      eff_U[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_U[i]->SetName(Form("eff_U_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_U[i]->SetDirectory(0);
    }
    recoEff_U_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_D_pass[i], recoEff_net[i]);
    recoEff_D[i] = (TH2D*)recoEff_D_pass[i]->Clone(Form("recoEff_D_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_D[i]->Divide(recoEff_net[i]);
    recoEff_D[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_D_pass[i]), *(recoEff_net[i]),"w") ){
      eff_D[i] = new TEfficiency(*(recoEff_D_pass[i]), *(recoEff_net[i]));
      eff_D[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_D[i]->SetName(Form("eff_D_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_D[i]->SetDirectory(0);
    }
    recoEff_D_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_photonU_pass[i], recoEff_net[i]);
    recoEff_photonU[i] = (TH2D*)recoEff_photonU_pass[i]->Clone(Form("recoEff_photonU_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_photonU[i]->Divide(recoEff_net[i]);
    recoEff_photonU[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_photonU_pass[i]), *(recoEff_net[i]),"w") ){
      eff_photonU[i] = new TEfficiency(*(recoEff_photonU_pass[i]), *(recoEff_net[i]));
      eff_photonU[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_photonU[i]->SetName(Form("eff_photonU_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_photonU[i]->SetDirectory(0);
    }

    recoEff_photonU_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_photonD_pass[i], recoEff_net[i]);
    recoEff_photonD[i] = (TH2D*)recoEff_photonD_pass[i]->Clone(Form("recoEff_photonD_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_photonD[i]->Divide(recoEff_net[i]);
    recoEff_photonD[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_photonD_pass[i]), *(recoEff_net[i]),"w") ){
      eff_photonD[i] = new TEfficiency(*(recoEff_photonD_pass[i]), *(recoEff_net[i]));
      eff_photonD[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_photonD[i]->SetName(Form("eff_photonD_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_photonD[i]->SetDirectory(0);
    }

    recoEff_photonD_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_noSF_pass[i], recoEff_noSF_net[i]);
    recoEff_noSF[i] = (TH2D*)recoEff_noSF_pass[i]->Clone(Form("recoEff_noSF_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_noSF[i]->Divide(recoEff_noSF_net[i]);
    recoEff_noSF[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_noSF_pass[i]), *(recoEff_noSF_net[i]),"w") ){
      eff_noSF[i] = new TEfficiency(*(recoEff_noSF_pass[i]), *(recoEff_noSF_net[i]));
      eff_noSF[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_noSF[i]->SetName(Form("eff_noSF_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_noSF[i]->SetDirectory(0);
    }
    recoEff_noSF_net[i]->SetDirectory(0);
    recoEff_noSF_pass[i]->SetDirectory(0);
  }
  
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_pt_pass[i], recoEff_pt_net[i]);
    recoEff_pt[i] = (TH1D*)recoEff_pt_pass[i]->Clone(Form("recoEff_pt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt[i]->Divide(recoEff_pt_net[i]);
    recoEff_pt[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pt_pass[i]), *(recoEff_pt_net[i]),"w") ){
      eff_pt[i] = new TEfficiency(*(recoEff_pt_pass[i]), *(recoEff_pt_net[i]));
      eff_pt[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_pt[i]->SetName(Form("eff_pt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_pt[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_pt_pass[i]->SetDirectory(0);
    recoEff_pt_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_pt_noPtWeight_pass[i], recoEff_pt_noPtWeight_net[i]);
    recoEff_pt_noPtWeight[i] = (TH1D*)recoEff_pt_noPtWeight_pass[i]->Clone(Form("recoEff_pt_noPtWeight_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_noPtWeight[i]->Divide(recoEff_pt_noPtWeight_net[i]);
    recoEff_pt_noPtWeight[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pt_noPtWeight_pass[i]), *(recoEff_pt_noPtWeight_net[i]),"w") ){
      eff_pt_noPtWeight[i] = new TEfficiency(*(recoEff_pt_noPtWeight_pass[i]), *(recoEff_pt_noPtWeight_net[i]));
      eff_pt_noPtWeight[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_pt_noPtWeight[i]->SetName(Form("eff_pt_noPtWeight_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_pt_noPtWeight[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_pt_noPtWeight_pass[i]->SetDirectory(0);
    recoEff_pt_noPtWeight_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_y_pass[i], recoEff_y_net[i]);
    recoEff_y[i] = (TH1D*)recoEff_y_pass[i]->Clone(Form("recoEff_y_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_y[i]->Divide(recoEff_y_net[i]);
    recoEff_y[i]->SetDirectory(0);
    
    forceConsistency(recoEff_yEE_pass[i], recoEff_yEE_net[i]);
    recoEff_yEE[i] = (TH1D*)recoEff_yEE_pass[i]->Clone(Form("recoEff_yEE_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_yEE[i]->Divide(recoEff_yEE_net[i]);
    recoEff_yEE[i]->SetDirectory(0);
    
    forceConsistency(recoEff_yBB_pass[i], recoEff_yBB_net[i]);
    recoEff_yBB[i] = (TH1D*)recoEff_yBB_pass[i]->Clone(Form("recoEff_yBB_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_yBB[i]->Divide(recoEff_yBB_net[i]);
    recoEff_yBB[i]->SetDirectory(0);
    
    forceConsistency(recoEff_yEB_pass[i], recoEff_yEB_net[i]);
    recoEff_yEB[i] = (TH1D*)recoEff_yEB_pass[i]->Clone(Form("recoEff_yEB_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_yEB[i]->Divide(recoEff_yEB_net[i]);
    recoEff_yEB[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_y_pass[i]), *(recoEff_y_net[i]),"w") ){
      eff_y[i] = new TEfficiency(*(recoEff_y_pass[i]), *(recoEff_y_net[i]));
      eff_y[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_y[i]->SetName(Form("eff_y_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_y[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_y_pass[i]->SetDirectory(0);
    recoEff_y_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_phi_pass[i], recoEff_phi_net[i]);
    recoEff_phi[i] = (TH1D*)recoEff_phi_pass[i]->Clone(Form("recoEff_phi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_phi[i]->Divide(recoEff_phi_net[i]);
    recoEff_phi[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_phi_pass[i]), *(recoEff_phi_net[i]),"w") ){
      eff_phi[i] = new TEfficiency(*(recoEff_phi_pass[i]), *(recoEff_phi_net[i]));
      eff_phi[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_phi[i]->SetName(Form("eff_phi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_phi[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_phi_pass[i]->SetDirectory(0);
    recoEff_phi_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_cent_pass[i], recoEff_cent_net[i]);
    recoEff_cent[i] = (TH1D*)recoEff_cent_pass[i]->Clone(Form("recoEff_cent_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_cent[i]->Divide(recoEff_cent_net[i]);
    recoEff_cent[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_cent_pass[i]), *(recoEff_cent_net[i]),"w") ){
      eff_cent[i] = new TEfficiency(*(recoEff_cent_pass[i]), *(recoEff_cent_net[i]));
      eff_cent[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_cent[i]->SetName(Form("eff_cent_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_cent[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_cent_pass[i]->SetDirectory(0);
    recoEff_cent_net[i]->SetDirectory(0);
  }

  forceConsistency(accept21_yields_pass[0], accept21_yields_net[0]);
  TEfficiency * accept_yields_TEff = new TEfficiency(*(accept21_yields_pass[0]), *(accept21_yields_net[0]));
  accept_yields_TEff->SetStatisticOption(TEfficiency::kBJeffrey);
  accept_yields_TEff->SetName("accept21_yields_TEff");
  accept_yields_TEff->SetDirectory(0);
  forceConsistency(accept21_y_pass[0], accept21_y_net[0]);
  TEfficiency * accept_y_TEff = new TEfficiency(*(accept21_y_pass[0]), *(accept21_y_net[0]));
  accept_y_TEff->SetStatisticOption(TEfficiency::kBJeffrey);
  accept_y_TEff->SetName("accept21_y_TEff");
  accept_y_TEff->SetDirectory(0);
  forceConsistency(accept21_pt_pass[0], accept21_pt_net[0]);
  TEfficiency * accept_pt_TEff = new TEfficiency(*(accept21_pt_pass[0]), *(accept21_pt_net[0]));
  accept_pt_TEff->SetStatisticOption(TEfficiency::kBJeffrey);
  accept_pt_TEff->SetName("accept21_pt_TEff");
  accept_pt_TEff->SetDirectory(0);

  zEff.~ZEfficiency();

  TFile * output;
  if(!isTest) output = new TFile(Form("resources/Z2ee_EfficiencyMC_%d.root",jobNumber),"recreate");
  else        output = new TFile(Form("resources/Z2ee_EfficiencyMC_%d_TEST.root",jobNumber),"recreate");
 
  for(int i = 0; i<nBins; i++){
    recoEff_net[i]->Write();
    recoEff_pass[i]->Write();
    recoEff[i]->Write();
    eff[i]->Write();
    
    recoEff_U_pass[i]->Write();
    recoEff_U[i]->Write();
    eff_U[i]->Write();
    
    recoEff_D_pass[i]->Write();
    recoEff_D[i]->Write();
    eff_D[i]->Write();
    
    recoEff_photonU[i]->Write();
    recoEff_photonU_pass[i]->Write();
    eff_photonU[i]->Write();
    
    recoEff_photonD[i]->Write();
    recoEff_photonD_pass[i]->Write();
    eff_photonD[i]->Write();
    
    recoEff_noSF_net[i]->Write();
    recoEff_noSF_pass[i]->Write();
    recoEff_noSF[i]->Write();
    eff_noSF[i]->Write();
    
    recoEff_pt[i]->Write();
    recoEff_pt_pass[i]->Write();
    recoEff_pt_net[i]->Write();
    eff_pt[i]->Write();
    
    recoEff_pt_noPtWeight[i]->Write();
    recoEff_pt_noPtWeight_pass[i]->Write();
    recoEff_pt_noPtWeight_net[i]->Write();
    eff_pt_noPtWeight[i]->Write();
    
    recoEff_y[i]->Write();
    recoEff_y_pass[i]->Write();
    recoEff_y_net[i]->Write();
    recoEff_y_net_smeared[i]->Write();
    eff_y[i]->Write();
    
    recoEff_yEE[i]->Write();
    recoEff_yEE_pass[i]->Write();
    recoEff_yEE_net[i]->Write();
    recoEff_yBB[i]->Write();
    recoEff_yBB_pass[i]->Write();
    recoEff_yBB_net[i]->Write();
    recoEff_yEB[i]->Write();
    recoEff_yEB_pass[i]->Write();
    recoEff_yEB_net[i]->Write();
    
    recoEff_phi[i]->Write();
    recoEff_phi_pass[i]->Write();
    recoEff_phi_net[i]->Write();
    eff_phi[i]->Write();
    
    recoEff_cent[i]->Write();
    recoEff_cent_pass[i]->Write();
    recoEff_cent_net[i]->Write();
    eff_cent[i]->Write();
    
    recoEff_pt_pass_forReso_Ratio_Reco[i] = (TH1D*) recoEff_pt_pass_forReso_Reco[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_Reco_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_Reco[i]->Divide(recoEff_pt_pass_forReso_Gen[i]);
    recoEff_pt_pass_forReso_Ratio_Reco[i]->Write();
    recoEff_pt_pass_forReso_Ratio_RecoSmeared[i] = (TH1D*) recoEff_pt_pass_forReso_RecoSmeared[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_RecoSmeared_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_RecoSmeared[i]->Divide(recoEff_pt_pass_forReso_Gen[i]);
    recoEff_pt_pass_forReso_Ratio_RecoSmeared[i]->Write();

    recoEff_pt_pass_forReso_Ratio_NominalToSmeared[i] = (TH1D*) recoEff_pt_pass_forReso_Ratio_RecoSmeared[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_SmearedToNominal_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_NominalToSmeared[i]->Divide(recoEff_pt_pass_forReso_Ratio_Reco[i]);
    recoEff_pt_pass_forReso_Ratio_NominalToSmeared[i]->Write();

    recoEff_pt_pass_forReso_Reco[i]->Write();
    recoEff_pt_pass_forReso_RecoSmeared[i]->Write();
    recoEff_pt_pass_forReso_Gen[i]->Write();
    recoEff_ptReso[i]->Write();

    massPeakOS[i]->Write();
    massPeakOS_noSF[i]->Write();
    massPeakSS[i]->Write();
    
    yReso[i]->Write();
  }
   
  yReco->Write();
  yGen->Write();
  yResponse->Write();
  yRatio = (TH1D*)yReco->Clone("yRatio");
  yRatio->Divide(yGen);
  yRatio->Write();
  
  TH1D * accept21_yields_ratio[120];
  TH1D * accept21_y_ratio[120];
  TH1D * accept21_pt_ratio[120];
  TH1D * accept21_scaleVariation_yields_ratio[120];
  TH1D * accept21_scaleVariation_y_ratio[120];
  TH1D * accept21_scaleVariation_pt_ratio[120];
  TH1D * accept21_nPDFVariation_yields_ratio[120];
  TH1D * accept21_nPDFVariation_y_ratio[120];
  TH1D * accept21_nPDFVariation_pt_ratio[120];
  TH1D * accept21_nPDFVariationMax_yields_ratio;
  TH1D * accept21_nPDFVariationMax_y_ratio;
  TH1D * accept21_nPDFVariationMax_pt_ratio;
  accept21_yields_noPtWeight_pass->Write();
  accept21_yields_noPtWeight_net->Write();
  accept21_pt_noPtWeight_pass->Write();
  accept21_pt_noPtWeight_net->Write();
  accept21_y_noPtWeight_pass->Write();
  accept21_y_noPtWeight_net->Write();
  for(int w = 0; w<weightHelper.getSize(); w++){
    accept21_yields_ratio[w] = (TH1D*)accept21_yields_pass[w]->Clone(Form("accept21_yields_ratio_%d",w));
    accept21_yields_ratio[w]->Divide(accept21_yields_net[w]);
    accept21_yields_ratio[w]->Write();
    
    accept21_y_ratio[w] = (TH1D*)accept21_y_pass[w]->Clone(Form("accept21_y_ratio_%d",w));
    accept21_y_ratio[w]->Divide(accept21_y_net[w]);
    accept21_y_ratio[w]->Write();
    
    accept21_pt_ratio[w] = (TH1D*)accept21_pt_pass[w]->Clone(Form("accept21_pt_ratio_%d",w));
    accept21_pt_ratio[w]->Divide(accept21_pt_net[w]);
    accept21_pt_ratio[w]->Write();

    //accept21_yields_pass[w]->Write();
    //accept21_yields_net[w]->Write();
    //accept21_pt_pass[w]->Write();
    //accept21_pt_net[w]->Write();
    //accept21_y_pass[w]->Write();
    //accept21_y_net[w]->Write();

    //scale variations
    if(w>=3 && w<=8){
      accept21_scaleVariation_yields_ratio[w] = (TH1D*) accept21_yields_ratio[w]->Clone(Form("accept21_scaleVariation_yields_%d",w));
      accept21_scaleVariation_yields_ratio[w]->Divide( accept21_yields_ratio[2] );
      accept21_scaleVariation_yields_ratio[w]->Write();
      accept21_scaleVariation_y_ratio[w] = (TH1D*) accept21_y_ratio[w]->Clone(Form("accept21_scaleVariation_y_%d",w));
      accept21_scaleVariation_y_ratio[w]->Divide( accept21_y_ratio[2] );
      accept21_scaleVariation_y_ratio[w]->Write();
      accept21_scaleVariation_pt_ratio[w] = (TH1D*) accept21_pt_ratio[w]->Clone(Form("accept21_scaleVariation_pt_%d",w));
      accept21_scaleVariation_pt_ratio[w]->Divide( accept21_pt_ratio[2] );
      accept21_scaleVariation_pt_ratio[w]->Write();
    }
    if(w==1 || w==2 || w>8){
      accept21_nPDFVariation_yields_ratio[w] = (TH1D*) accept21_yields_ratio[w]->Clone(Form("accept21_nPDFVariation_yields_%d",w));
      accept21_nPDFVariation_yields_ratio[w]->Divide( accept21_yields_ratio[0] );
      accept21_nPDFVariation_yields_ratio[w]->Write();
      accept21_nPDFVariation_y_ratio[w] = (TH1D*) accept21_y_ratio[w]->Clone(Form("accept21_nPDFVariation_y_%d",w));
      accept21_nPDFVariation_y_ratio[w]->Divide( accept21_y_ratio[0] );
      accept21_nPDFVariation_y_ratio[w]->Write();
      accept21_nPDFVariation_pt_ratio[w] = (TH1D*) accept21_pt_ratio[w]->Clone(Form("accept21_nPDFVariation_pt_%d",w));
      accept21_nPDFVariation_pt_ratio[w]->Divide( accept21_pt_ratio[0] );
      accept21_nPDFVariation_pt_ratio[w]->Write();
 
      if(w>8){
        if(w==9){
          accept21_nPDFVariationMax_yields_ratio = (TH1D*) accept21_nPDFVariation_yields_ratio[w]->Clone("accept21_nPDFVariationMax_yields_ratio");
          accept21_nPDFVariationMax_yields_ratio->Reset();
          accept21_nPDFVariationMax_y_ratio = (TH1D*) accept21_nPDFVariation_y_ratio[w]->Clone("accept21_nPDFVariationMax_y_ratio");
          accept21_nPDFVariationMax_y_ratio->Reset();
          accept21_nPDFVariationMax_pt_ratio = (TH1D*) accept21_nPDFVariation_pt_ratio[w]->Clone("accept21_nPDFVariationMax_pt_ratio");
          accept21_nPDFVariationMax_pt_ratio->Reset();
        }
        for(int bin = 0; bin<accept21_nPDFVariation_yields_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept21_nPDFVariation_yields_ratio[w]->GetBinContent(bin) );
          if( content > accept21_nPDFVariationMax_yields_ratio->GetBinContent(bin)){
            accept21_nPDFVariationMax_yields_ratio->SetBinContent(bin, content);
          }
        }
        for(int bin = 0; bin<accept21_nPDFVariation_y_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept21_nPDFVariation_y_ratio[w]->GetBinContent(bin) );
          if( content > accept21_nPDFVariationMax_y_ratio->GetBinContent(bin)){
            accept21_nPDFVariationMax_y_ratio->SetBinContent(bin, content);
          }
        }
        for(int bin = 0; bin<accept21_nPDFVariation_pt_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept21_nPDFVariation_pt_ratio[w]->GetBinContent(bin) );
          if( content > accept21_nPDFVariationMax_pt_ratio->GetBinContent(bin)){
            accept21_nPDFVariationMax_pt_ratio->SetBinContent(bin, content);
          }
        }
      }
    }
  }
  accept21_nPDFVariationMax_yields_ratio->Scale(1.0/1.645);//envelope gives 90%, scale down to 68% coverage
  accept21_nPDFVariationMax_y_ratio->Scale(1.0/1.645);
  accept21_nPDFVariationMax_pt_ratio->Scale(1.0/1.645);
  accept21_nPDFVariationMax_yields_ratio->Write();
  accept21_nPDFVariationMax_y_ratio->Write();
  accept21_nPDFVariationMax_pt_ratio->Write();
  
  TH1D * accept21_pdfTypeVariation_yields_ratio = (TH1D*) accept21_nPDFVariation_yields_ratio[1]->Clone("accept21_pdfTypeVariation_yields_ratio");
  TH1D * accept21_totalUncert_yields = (TH1D*) accept21_nPDFVariationMax_yields_ratio->Clone("accept21_totalUncert_yields");
  for(int bin = 0; bin<accept21_pdfTypeVariation_yields_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept21_pdfTypeVariation_yields_ratio->GetBinContent(bin) );
    accept21_pdfTypeVariation_yields_ratio->SetBinContent(bin, content);

    float content2 = accept21_totalUncert_yields->GetBinContent(bin);
    accept21_totalUncert_yields->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept21_pdfTypeVariation_yields_ratio->Write();
  accept21_totalUncert_yields->Write(); 
 
  TH1D * accept21_pdfTypeVariation_y_ratio = (TH1D*) accept21_nPDFVariation_y_ratio[1]->Clone("accept21_pdfTypeVariation_y_ratio");
  TH1D * accept21_totalUncert_y = (TH1D*) accept21_nPDFVariationMax_y_ratio->Clone("accept21_totalUncert_y");
  for(int bin = 0; bin<accept21_pdfTypeVariation_y_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept21_pdfTypeVariation_y_ratio->GetBinContent(bin) );
    accept21_pdfTypeVariation_y_ratio->SetBinContent(bin, content);
    
    float content2 = accept21_totalUncert_y->GetBinContent(bin);
    accept21_totalUncert_y->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept21_pdfTypeVariation_y_ratio->Write();
  accept21_totalUncert_y->Write();  

  TH1D * accept21_pdfTypeVariation_pt_ratio = (TH1D*) accept21_nPDFVariation_pt_ratio[1]->Clone("accept21_pdfTypeVariation_pt_ratio");
  TH1D * accept21_totalUncert_pt = (TH1D*) accept21_nPDFVariationMax_pt_ratio->Clone("accept21_totalUncert_pt");
  for(int bin = 0; bin<accept21_pdfTypeVariation_pt_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept21_pdfTypeVariation_pt_ratio->GetBinContent(bin) );
    accept21_pdfTypeVariation_pt_ratio->SetBinContent(bin, content);
    
    float content2 = accept21_totalUncert_pt->GetBinContent(bin);
    accept21_totalUncert_pt->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept21_pdfTypeVariation_pt_ratio->Write();
  accept21_totalUncert_pt->Write(); 

  accept_yields_TEff->Write();
  accept_y_TEff->Write();
  accept_pt_TEff->Write();
  
  unfolding_recoPt->Write();
  unfolding_genPt->Write();
  unfolding_response->Write();
  unfolding_response_noPtWeight->Write();
  unfolding_recoPtU->Write();
  unfolding_genPtU->Write();
  unfolding_responseU->Write();
  unfolding_recoPtD->Write();
  unfolding_genPtD->Write();
  unfolding_responseD->Write();

  float Asq = 208*208;
  //float TAA090 = 6.274;
  //float TAA090RelError = 0.022;
  //float TAA090RelError = 0.0;//do not include TAA for now
  float unitConversion = TMath::Power(10,9);
  for(int w = 282; w<=338; w++){
    theoryPtCT14[w-282]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryPtCT14[w-282] );
  }
  for(int w = 282; w<=338; w++){
    theoryRapCT14[w-282]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryRapCT14[w-282] );
  }
  for(int w = 1177; w<=1209; w++){
    theoryPtnCTEQ15[w-1177]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryPtnCTEQ15[w-1177] );
  }
  for(int w = 1177; w<=1209; w++){
    theoryRapnCTEQ15[w-1177]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryRapnCTEQ15[w-1177] );
  }
  for(int w = 1080; w<=1176; w++){
    theoryPtEPPS16[w-1080]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryPtEPPS16[w-1080] );
  }
  for(int w = 1080; w<=1176; w++){
    theoryRapEPPS16[w-1080]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryRapEPPS16[w-1080] );
  }
  for(int w = 1080; w<=1176; w++){
    theoryNetEPPS16[w-1080]->Scale(s.DY_XS*0.5/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
  }
  for(int w = 0; w<=8; w++){
    theoryPtScaleVariations[w]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryPtScaleVariations[w] );
    theoryPtScaleVariations[w]->Write();
    if(w>0) theoryPtScaleVariations[w]->Divide(theoryPtScaleVariations[0]);
  }
  for(int w = 0; w<=8; w++){
    theoryRapScaleVariations[w]->Scale(s.DY_XS*0.5*Asq/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    histHelper.makeDifferential( theoryRapScaleVariations[w] );
    theoryRapScaleVariations[w]->Write();
    if(w>0) theoryRapScaleVariations[w]->Divide(theoryRapScaleVariations[0]);
  }
  for(int w = 0; w<=8; w++){
    theoryNetScaleVariations[w]->Scale(s.DY_XS*0.5/theoryEventWeights.at(0)/unitConversion); //add the factor of 0.5 because we use both the mu and e channel for this, but for comparison we only use one channel
    theoryNetScaleVariations[w]->Write();
    if(w>0) theoryNetScaleVariations[w]->Divide(theoryNetScaleVariations[0]);
  }

  TH1D* theoryPtCT14_Band = new TH1D("theoryCT14_pt_Band","theoryCT14_pt_Band",s.nZPtBins-1,s.zPtBins);
  TH1D* theoryPtnCTEQ15_Band = new TH1D("theorynCTEQ15_pt_Band","theorynCTEQ15_pt_Band",s.nZPtBins-1,s.zPtBins);
  TH1D* theoryPtEPPS16_Band = new TH1D("theoryEPPS16_pt_Band","theoryEPPS16_pt_Band",s.nZPtBins-1,s.zPtBins);
  TH1D* theoryRapCT14_Band = new TH1D("theoryCT14_rap_Band","theoryCT14_rap_Band",s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D* theoryRapnCTEQ15_Band = new TH1D("theorynCTEQ15_rap_Band","theorynCTEQ15_rap_Band",s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D* theoryRapEPPS16_Band = new TH1D("theoryEPPS16_rap_Band","theoryEPPS16_rap_Band",s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D* theoryNetEPPS16_Band = new TH1D("theoryEPPS16_Net_Band","theoryEPPS16_Net_Band",1,0,1);
  TH1D* theoryRapScales_Band = new TH1D("theoryScales_rap_Band","theoryScales_rap_Band",s.nZRapBins,-s.maxZRap,s.maxZRap);
  TH1D* theoryNetScales_Band = new TH1D("theoryScales_Net_Band","theoryScales_Net_Band",1,0,1);
  TH1D* theoryPtScales_Band = new TH1D("theoryScales_pt_Band","theoryScales_pt_Band",s.nZPtBins-1,s.zPtBins);
  float max[100];
  float min[100];
  //Scales
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1; w<=8; w++){
    for(int i = 0; i<theoryPtScales_Band->GetSize(); i++){
      if(theoryPtScaleVariations[w]->GetBinContent(i) < min[i]) min[i] = theoryPtScaleVariations[w]->GetBinContent(i);
      if(theoryNetScaleVariations[w]->GetBinContent(i) > max[i]) max[i] = theoryPtScaleVariations[w]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryPtScales_Band->GetSize(); i++){
    theoryPtScales_Band->SetBinContent(i,1);
    theoryPtScales_Band->SetBinError(i,TMath::Max(max[i]-1,1-min[i]) );//scales envelope symmetrized
  }
  theoryPtScales_Band->Write();
  
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1; w<=8; w++){
    for(int i = 0; i<theoryRapScales_Band->GetSize(); i++){
      if(theoryRapScaleVariations[w]->GetBinContent(i) < min[i]) min[i] = theoryRapScaleVariations[w]->GetBinContent(i);
      if(theoryNetScaleVariations[w]->GetBinContent(i) > max[i]) max[i] = theoryRapScaleVariations[w]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryRapScales_Band->GetSize(); i++){
    theoryRapScales_Band->SetBinContent(i,1);
    theoryRapScales_Band->SetBinError(i,TMath::Max(max[i]-1,1-min[i]) );//scales envelope symmetrized
  }
  theoryRapScales_Band->Write();
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1; w<=8; w++){
    for(int i = 0; i<theoryNetScales_Band->GetSize(); i++){
      if(theoryNetScaleVariations[w]->GetBinContent(i) < min[i]) min[i] = theoryNetScaleVariations[w]->GetBinContent(i);
      if(theoryNetScaleVariations[w]->GetBinContent(i) > max[i]) max[i] = theoryNetScaleVariations[w]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryNetScales_Band->GetSize(); i++){
    theoryNetScales_Band->SetBinContent(i,1);
    theoryNetScales_Band->SetBinError(i,TMath::Max(max[i]-1,1-min[i]) );//scales envelope symmetrized
  }
  theoryNetScales_Band->Write();

  //CT14
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 282; w<=338; w++){
    for(int i = 0; i<theoryPtCT14_Band->GetSize(); i++){
      if(theoryPtCT14[w-282]->GetBinContent(i) < min[i]) min[i] = theoryPtCT14[w-282]->GetBinContent(i);
      if(theoryPtCT14[w-282]->GetBinContent(i) > max[i]) max[i] = theoryPtCT14[w-282]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryPtCT14_Band->GetSize(); i++){
    theoryPtCT14_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryPtCT14_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2) + theoryPtScales_Band->GetBinError(i)*theoryPtScales_Band->GetBinError(i)*theoryPtCT14_Band->GetBinContent(i)*theoryPtCT14_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryPtCT14_Band->Write();
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 282; w<=338; w++){
    for(int i = 0; i<theoryRapCT14_Band->GetSize(); i++){
      if(theoryRapCT14[w-282]->GetBinContent(i) < min[i]) min[i] = theoryRapCT14[w-282]->GetBinContent(i);
      if(theoryRapCT14[w-282]->GetBinContent(i) > max[i]) max[i] = theoryRapCT14[w-282]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryRapCT14_Band->GetSize(); i++){
    theoryRapCT14_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryRapCT14_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2)  + theoryRapScales_Band->GetBinError(i)*theoryRapScales_Band->GetBinError(i)*theoryRapCT14_Band->GetBinContent(i)*theoryRapCT14_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryRapCT14_Band->Write();
 
  //nCTEQ15 
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1177; w<=1209; w++){
    for(int i = 0; i<theoryPtnCTEQ15_Band->GetSize(); i++){
      if(theoryPtnCTEQ15[w-1177]->GetBinContent(i) < min[i]) min[i] = theoryPtnCTEQ15[w-1177]->GetBinContent(i);
      if(theoryPtnCTEQ15[w-1177]->GetBinContent(i) > max[i]) max[i] = theoryPtnCTEQ15[w-1177]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryPtnCTEQ15_Band->GetSize(); i++){
    theoryPtnCTEQ15_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryPtnCTEQ15_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2) + theoryPtScales_Band->GetBinError(i)*theoryPtScales_Band->GetBinError(i)*theoryPtnCTEQ15_Band->GetBinContent(i)*theoryPtnCTEQ15_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryPtnCTEQ15_Band->Write();
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1177; w<=1209; w++){
    for(int i = 0; i<theoryRapnCTEQ15_Band->GetSize(); i++){
      if(theoryRapnCTEQ15[w-1177]->GetBinContent(i) < min[i]) min[i] = theoryRapnCTEQ15[w-1177]->GetBinContent(i);
      if(theoryRapnCTEQ15[w-1177]->GetBinContent(i) > max[i]) max[i] = theoryRapnCTEQ15[w-1177]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryRapnCTEQ15_Band->GetSize(); i++){
    theoryRapnCTEQ15_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryRapnCTEQ15_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2)  + theoryRapScales_Band->GetBinError(i)*theoryRapScales_Band->GetBinError(i)*theoryRapnCTEQ15_Band->GetBinContent(i)*theoryRapnCTEQ15_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryRapnCTEQ15_Band->Write();

  //EPPS16
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1080; w<=1176; w++){
    for(int i = 0; i<theoryPtEPPS16_Band->GetSize(); i++){
      if(theoryPtEPPS16[w-1080]->GetBinContent(i) < min[i]) min[i] = theoryPtEPPS16[w-1080]->GetBinContent(i);
      if(theoryPtEPPS16[w-1080]->GetBinContent(i) > max[i]) max[i] = theoryPtEPPS16[w-1080]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryPtEPPS16_Band->GetSize(); i++){
    theoryPtEPPS16_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryPtEPPS16_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2) + theoryPtScales_Band->GetBinError(i)*theoryPtScales_Band->GetBinError(i)*theoryPtEPPS16_Band->GetBinContent(i)*theoryPtEPPS16_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryPtEPPS16_Band->Write();
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1080; w<=1176; w++){
    for(int i = 0; i<theoryRapEPPS16_Band->GetSize(); i++){
      if(theoryRapEPPS16[w-1080]->GetBinContent(i) < min[i]) min[i] = theoryRapEPPS16[w-1080]->GetBinContent(i);
      if(theoryRapEPPS16[w-1080]->GetBinContent(i) > max[i]) max[i] = theoryRapEPPS16[w-1080]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryRapEPPS16_Band->GetSize(); i++){
    theoryRapEPPS16_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryRapEPPS16_Band->SetBinError(i,TMath::Sqrt( TMath::Power( ((max[i]-min[i])/2.0/1.644854),2)  + theoryRapScales_Band->GetBinError(i)*theoryRapScales_Band->GetBinError(i)*theoryRapEPPS16_Band->GetBinContent(i)*theoryRapEPPS16_Band->GetBinContent(i) ) );//scaled down to 1 sigma level and add TAA
  }
  theoryRapEPPS16_Band->Write();
  for(int i = 0; i<100; i++){ max[i] = 0; min[i] = 9999999;}
  for(int w = 1080; w<=1176; w++){
    for(int i = 0; i<theoryNetEPPS16_Band->GetSize(); i++){
      if(theoryNetEPPS16[w-1080]->GetBinContent(i) < min[i]) min[i] = theoryNetEPPS16[w-1080]->GetBinContent(i);
      if(theoryNetEPPS16[w-1080]->GetBinContent(i) > max[i]) max[i] = theoryNetEPPS16[w-1080]->GetBinContent(i);
    }
  }
  for(int i = 0; i<theoryNetEPPS16_Band->GetSize(); i++){
    theoryNetEPPS16_Band->SetBinContent(i,((max[i]+min[i])/2.0));
    theoryNetEPPS16_Band->SetBinError(i,TMath::Sqrt(TMath::Power((max[i]-min[i])/2.0/1.644854,2)  + theoryNetScales_Band->GetBinError(i)*theoryNetScales_Band->GetBinError(i)*theoryNetEPPS16_Band->GetBinContent(i)*theoryNetEPPS16_Band->GetBinContent(i) )) ;//scaled down to 1 sigma level and add TAA
  }
  theoryNetEPPS16_Band->Write();

  pdfXVsPt->Write();
  pdfXVsRap->Write();

  for(int i = 0; i<12; i++){
    pdfIDVsPt[i]->Divide(pdfIDVsPt[12]);
    pdfIDVsPt[i]->Write();
  }

  pdfX1X2LowPt->Write();
  pdfX1X2HighPt->Write();

  
  float diagonals[8] = {0};
  for(int i = 0; i<8; i++){
    diagonals[i] = TMath::Sqrt(covariance->GetBinContent(i+1,i+1));
  }  
  for(int i = 0; i<8; i++){
    for(int j = 0; j<8; j++){
      covariance->SetBinContent(i+1,j+1, covariance->GetBinContent(i+1,j+1)/diagonals[i]/diagonals[j]);
    }
  }
  covariance->Write();

  output->Close();

  timer.Stop();
  timer.Report();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: Z_EE_EfficiencyMC <fileList> <isTest> <job #> <total number of jobs>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  bool isTest = (bool)std::atoi(argv[2]);
  int job = std::atoi(argv[3]);
  int totalJobs = std::atoi(argv[4]);

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  doZ2EE(listOfFiles, job, isTest);
  return 0; 
}


int pdfProcess(int id1,int id2,int nPartons){
  if(nPartons>2) return 11;//NNNLO or higher

  //uubar
  if((id1==2 && id2==-2) || (id1==-2 && id2==2)){
    if(nPartons==0) return 0;
    if(nPartons==1) return 3;
  }
  //ddbar
  if((id1==1 && id2==-1) || (id1==-1 && id2==1)){
    if(nPartons==0) return 1;
    if(nPartons==1) return 4;
  }
  //ssbar + QQbar
  if((id1==3 && id2==-3) || (id1==-3 && id2==3) || (id1==4 && id2==-4) || (id1==-4 && id2==4) || (id1==5 && id2==-5) || (id1==-5 && id2==5) || (id1==6 && id2==-6) || (id1==-6 && id2==6)){
    if(nPartons==0) return 2;
    if(nPartons==1) return 5;
  }
  //gu
  if((id1==21 && id2==2) || (id1==2 && id2==21)){
    if(nPartons==1) return 6;
  }
  //gd
  if((id1==21 && id2==1) || (id1==1 && id2==21)){
    if(nPartons==1) return 7;
  }
  //gQ
  if((id1==21) || (id2==21)){
    if(nPartons==1 || nPartons==0) return 8;//NLO
    if(nPartons==2) return 9;//NNLO
  }else{
    if(nPartons==2) return 10;//NNLO no gluons initial state 
  }
  return 11;
}

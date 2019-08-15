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
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2EE(std::vector< std::string > files, int jobNumber, bool isMC, std::string outputTag,int hiBinVar = 0, bool doZDC = false){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");


  HistNameHelper h = HistNameHelper();
  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronEnergyScale energyScale = ElectronEnergyScale("data");
  ElectronEnergyScale energyScaleMC = ElectronEnergyScale("MC");
  Settings s = Settings();
  
  ZEfficiency zEff = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", isMC);
  ZEfficiency zEffU = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", isMC, 1);
  ZEfficiency zEffD = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", isMC, -1);
  ZEfficiency zEffphotonU = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", isMC, 2);
  ZEfficiency zEffphotonD = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root", isMC, -2);
  
  MCReweight * vzRW;
  if(isMC) vzRW = new MCReweight("resources/vzReweight.root");

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();
  
  TH1D * nEvents = new TH1D("nEvents","nEvents",1,0,1);

  TH2D * dEtaVsdEtaSeed = new TH2D("dEtaVsdEtaSeed","dEtaVsdEtaSeed",50,0,0.01,50,0,0.01);

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TH1D * massPeakOS_withEff[nBins]; 
  TH1D * massPeakSS_withEff[nBins]; 
  TH2D * massVsPt[nBins];
  TH1D * yields;
  TH1D * yieldsSS;
  
  TH1D * massPeakOS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * massPeakSS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * massPeakOS_ptLT0p5_withEff[nBins]; 
  TH1D * massPeakSS_ptLT0p5_withEff[nBins]; 
  TH1D * massPeakOS_TauTau_withEff[nBins]; 
  TH1D * massPeakSS_TauTau_withEff[nBins]; 
  
  TH1D * pTOS_withEff[nBins][5]; 
  TH1D * pTSS_withEff[nBins][5]; 
  TH1D * pTOS_ptLT0p5acoLT0p001_withEff[nBins][5]; 
  TH1D * pTOS_TauTau_withEff[nBins]; 
  TH1D * pTSS_TauTau_withEff[nBins]; 
  
  TH1D * yOS_withEff[nBins][5]; 
  TH1D * ySS_withEff[nBins][5]; 
  TH1D * yOS_ptLT0p5acoLT0p001_withEff[nBins][5]; 
  TH1D * yOS_TauTau_withEff[nBins]; 
  TH1D * ySS_TauTau_withEff[nBins]; 
  
  TH1D * yieldOS_withEff[nBins][5]; 
  TH1D * yieldSS_withEff[nBins][5]; 
  TH1D * yieldOS_ptLT0p5acoLT0p001_withEff[nBins][5]; 
  TH1D * yieldOS_TauTau_withEff[nBins]; 
  TH1D * yieldSS_TauTau_withEff[nBins]; 
  
  TProfile * v2Num[nBins];
  TProfile * v2NumVsCent[3];
  TProfile * v2Denom[nBins];
  TProfile * v2DenomVsCent;
  TProfile * v2Q1Mid[nBins];
  TProfile * v2Q1MidVsCent;
  TProfile * v2Q2Mid[nBins];
  TProfile * v2Q2MidVsCent;
  TProfile * v2AvgEffVsCent[3];
  
  TProfile * v2EleNum[nBins];
  TProfile * v2EleNumVsCent;
  TProfile * v2EleDenom[nBins];
  TProfile * v2EleDenomVsCent;
  TProfile * v2EleQ1Mid[nBins];
  TProfile * v2EleQ1MidVsCent;
  TProfile * v2EleQ2Mid[nBins];
  TProfile * v2EleQ2MidVsCent;

  TH1D * candPt[nBins];
  TH1D * candPtFine[nBins];
  TH1D * candPtFiner[nBins];
  TH1D * candPt_unnormalized[nBins];
  TH1D * candEta[nBins];
  TH1D * candY[nBins];
  TH1D * candPhi[nBins];
  TH2D * candPtVsM[nBins];
  TH2D * candAcoVsM[nBins];
  TH2D * candAcoVsPt[nBins];
  TH1D * candAco[nBins];
  
  TH1D * lepPt;  
  TH1D * lepEta;
  TH1D * lepPhi;  


  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_withEff[i] = new TH1D(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_withEff[i] = new TH1D(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massVsPt[i] = new TH2D(Form("massVsPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}};p_{T}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1],s.nZPtBins-1,s.zPtBins);

    massPeakOS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("massPeakOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("massPeakSS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_ptLT0p5_withEff[i] = new TH1D(Form("massPeakOS_ptLT0p5_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_ptLT0p5_withEff[i] = new TH1D(Form("massPeakSS_ptLT0p5_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_TauTau_withEff[i] = new TH1D(Form("massPeakOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_TauTau_withEff[i] = new TH1D(Form("massPeakSS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    for(int j = 0; j < (isMC ? 1 : 5); j++){
      pTOS_withEff[i][j] = new TH1D(Form("pTOS_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
      pTSS_withEff[i][j] = new TH1D(Form("pTSS_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
      pTOS_ptLT0p5acoLT0p001_withEff[i][j] = new TH1D(Form("pTOS_ptLT0p5acoLT0p001_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    }
    pTOS_TauTau_withEff[i] = new TH1D(Form("pTOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    pTSS_TauTau_withEff[i] = new TH1D(Form("pTSS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);

    int nRapBins = 14;
    for(int j = 0; j < (isMC ? 1 : 5); j++){
      yOS_withEff[i][j] = new TH1D(Form("yOS_withEff%s_%d_%d",h.variationName.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-s.maxZRapEle,s.maxZRapEle);
      ySS_withEff[i][j] = new TH1D(Form("ySS_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-s.maxZRapEle,s.maxZRapEle);
      yOS_ptLT0p5acoLT0p001_withEff[i][j] = new TH1D(Form("yOS_ptLT0p5acoLT0p001_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-s.maxZRapEle,s.maxZRapEle);
    }
    yOS_TauTau_withEff[i] = new TH1D(Form("yOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-s.maxZRapEle,s.maxZRapEle);
    ySS_TauTau_withEff[i] = new TH1D(Form("ySS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-s.maxZRapEle,s.maxZRapEle);
    
    for(int j = 0; j < (isMC ? 1 : 5); j++){
      yieldOS_withEff[i][j] = new TH1D(Form("yieldOS_withEff%s_%d_%d",h.variationName.at(j).c_str(),c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
      yieldSS_withEff[i][j] = new TH1D(Form("yieldSS_withEff%s_%d_%d",h.variationName.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",1,0,1);
      yieldOS_ptLT0p5acoLT0p001_withEff[i][j] = new TH1D(Form("yieldOS_ptLT0p5acoLT0p001_withEff%s_%d_%d",h.variationName.at(j).c_str(), c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
    }
    yieldOS_TauTau_withEff[i] = new TH1D(Form("yieldOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
    yieldSS_TauTau_withEff[i] = new TH1D(Form("yieldSS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);

    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q1Mid[i] = new TProfile(Form("v2Q1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q2Mid[i] = new TProfile(Form("v2Q2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
    v2EleNum[i] = new TProfile(Form("v2EleNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleDenom[i] = new TProfile(Form("v2EleDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleQ1Mid[i] = new TProfile(Form("v2EleQ1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleQ2Mid[i] = new TProfile(Form("v2EleQ2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
    candPt[i] = new TH1D(Form("candPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candPtFine[i] = new TH1D(Form("candPtFine_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",400,0,200);
    candPtFiner[i] = new TH1D(Form("candPtFiner_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",800,0,200);
    candPt_unnormalized[i] = new TH1D(Form("candPt_unnormalized_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candEta[i] = new TH1D(Form("candEta_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",14,-s.maxZRapEle,s.maxZRapEle);
    candY[i] = new TH1D(Form("candY_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",14,-s.maxZRapEle,s.maxZRapEle);
    candPhi[i] = new TH1D(Form("candPhi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-TMath::Pi(),TMath::Pi());
    candPtVsM[i] = new TH2D(Form("candPtVsM_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",40,0,10,60,60,120);
    candAcoVsM[i] = new TH2D(Form("candAcoVsM_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",50,0,0.005,60,60,120); 
    candAcoVsPt[i] = new TH2D(Form("candAcoVsPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",50,0,0.005,50,0,5);
    candAco[i] = new TH1D(Form("candAco_%d_%d",c.getCentBinLow(i), c.getCentBinHigh(i)),"",50,0,0.005); 
  }  
  yields = new TH1D("yields","yields",nBins,0,nBins);  
  yieldsSS = new TH1D("yieldsSS","yieldsSS",nBins,0,nBins);  
  
  lepPt = new TH1D("lepPt",";p_{T}",s.nZPtBins-1,s.zPtBins);
  lepEta = new TH1D("lepEta",";p_{T}",20,-s.maxZRap,s.maxZRap);
  lepPhi = new TH1D("lepPhi",";p_{T}",20,-TMath::Pi(),TMath::Pi());

  for(int i = 0; i<3; i++){ 
    v2AvgEffVsCent[i] = new TProfile(Form("v2AvgEffVsCent%s",h.variationName.at(i).c_str()),"",nBins,0,nBins);
    v2NumVsCent[i] = new TProfile(Form("v2NumVsCent%s",h.variationName.at(i).c_str()),"",nBins,0,nBins);
  }
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2Q1MidVsCent = new TProfile("v2Q1MidVsCent","",nBins,0,nBins);
  v2Q2MidVsCent = new TProfile("v2Q2MidVsCent","",nBins,0,nBins);

  v2EleNumVsCent = new TProfile("v2EleNumVsCent","",nBins,0,nBins);
  v2EleDenomVsCent = new TProfile("v2EleDenomVsCent","",nBins,0,nBins);
  v2EleQ1MidVsCent = new TProfile("v2EleQ1MidVsCent","",nBins,0,nBins);
  v2EleQ2MidVsCent = new TProfile("v2EleQ2MidVsCent","",nBins,0,nBins);

  int nEle;
  int hiBin;
  int hiBinZDC;
  int hiNpix;
  float hiHF;
  float vz;

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
  std::vector< float > * eledEtaAtVtx = 0;
  std::vector< float > * eledPhiAtVtx = 0;
  std::vector< float > * eleSCEta = 0;
  std::vector< float > * eleSeedEta = 0;
  std::vector< float > * eleSCPhi = 0;
  std::vector< float > * eleHoverEBc = 0;
  //std::vector< float > * eleD0 = 0;
  //std::vector< float > * eleDz = 0;
  std::vector< float > * eleIP3D = 0;
  std::vector< float > * eleEoverPInv = 0;

  int nMC;
  std::vector< int > * mcPID = 0;
  std::vector< int > * mcMomPID = 0;
  std::vector< int > * mcGMomPID = 0;
  

  int ZDC_n = 0;
  float ZDC_E[100];

  int hiNevtPlane;
  float hiQVecMag[200];
  float hiQVecAngle[200];

  for(unsigned int f = 0; f<files.size(); f++){
    timer.StartSplit("Opening Files");

    TFile * in = TFile::Open(files.at(f).c_str(),"read");
    if(f%5 == 0)  std::cout << f << "/" << files.size() << std::endl;

    TTree * hltTree = (TTree*)in->Get("hltanalysis/HltTree");
    hltTree->SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v1",&HLT_DoubleEle10); 
    hltTree->SetBranchAddress("HLT_HIEle20Gsf_v1",&HLT_SingleEle20); 

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    eTree->SetBranchAddress("nEle",&nEle);
    eTree->SetBranchAddress("elePt",&elePt);
    eTree->SetBranchAddress("eleEta",&eleEta);
    eTree->SetBranchAddress("elePhi",&elePhi);
    eTree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta);
    eTree->SetBranchAddress("eleMissHits",&eleMissHits);
    eTree->SetBranchAddress("eleCharge",&eleCharge);
    eTree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
    eTree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
    eTree->SetBranchAddress("eleSCEta",&eleSCEta);
    eTree->SetBranchAddress("eleSeedEta",&eleSeedEta);
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverEBc",&eleHoverEBc);
    //eTree->SetBranchAddress("eleD0",&eleD0);
    //eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleIP3D",&eleIP3D);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);
    
    if(isMC){
      eTree->SetBranchAddress("nMC",&nMC);
      eTree->SetBranchAddress("mcPID",&mcPID);
      eTree->SetBranchAddress("mcMomPID",&mcMomPID);
      eTree->SetBranchAddress("mcGMomPID",&mcGMomPID);
    }

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    //evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("hiNpix",&hiNpix);
    evtTree->SetBranchAddress("hiHF",&hiHF);
    evtTree->SetBranchAddress("vz",&vz);
    evtTree->SetBranchAddress("hiNevtPlane",&hiNevtPlane);  
    evtTree->SetBranchAddress("hiQVecMag",hiQVecMag);  
    evtTree->SetBranchAddress("hiQVecAngle",hiQVecAngle);  
    if(isMC) evtTree->SetBranchAddress("ttbar_w",ttbar_w);
  
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

    //TTree * L1Tree = (TTree*)in->Get("l1object/L1UpgradeFlatTree");
    //L1Tree->SetBranchAddress("nEGs",&(eTrig.L1nEGs));
    //L1Tree->SetBranchAddress("egEta", &(eTrig.L1egEta));
    //L1Tree->SetBranchAddress("egPhi", &(eTrig.L1egPhi));
    //L1Tree->SetBranchAddress("egEt", &(eTrig.L1egEt));

    TTree * ZDCTree;
    ZDCTree = (TTree*)in->Get("rechitanalyzerpp/zdcrechit");
    ZDCTree->SetBranchAddress("n",&ZDC_n);
    ZDCTree->SetBranchAddress("e",ZDC_E);

    TTree * HLTObjTree;
    HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIEle20Gsf_v");
    HLTObjTree->SetBranchAddress("eta",&(eTrig.HLTEta));
    HLTObjTree->SetBranchAddress("phi",&(eTrig.HLTPhi));
    HLTObjTree->SetBranchAddress("pt",&(eTrig.HLTPt));

    bool areAllEleBranchesOn = true;
    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      timer.StartSplit("Checking Evt Selections");
      skimTree->GetEntry(i);
      if(! (pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter)) continue;
      
      //event selections
      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;
      hiBin = cb.getHiBinFromhiHF(hiHF, (isMC ? 3 : hiBinVar) );
    
      if(doZDC) ZDCTree->GetEntry(i);
      float ZDCSum = 0;
      for(int z = 0; z< ZDC_n; z++) ZDCSum += ZDC_E[z];
      if(!isMC && doZDC) hiBinZDC = cb.getHiBinFromZDC( hiNpix , ZDCSum , 0);
      else hiBinZDC = hiBin;
      

      float eventWeight = 1.0;
      if(isMC) eventWeight = vzRW->reweightFactor( vz ) * c.findNcoll( hiBin ) * (ttbar_w->at(1080)/10000.0);//1080 is EEPS16NLO+CT14
      nEvents->Fill(0.5,eventWeight);
      
      timer.StartSplit("Checking Number of electrons");
      //check if the event has at least 2 electrons
      //some clever branch status changing is done here to avoid loading too much stuff for the simple check
      if(areAllEleBranchesOn){
        eTree->SetBranchStatus("*",0);
        eTree->SetBranchStatus("nEle",1);
        areAllEleBranchesOn = false;
      }
      eTree->GetEntry(i);
      if(nEle<2) continue;
      
      timer.StartSplit("Checking HLT Selections");
      //check for the trigger we want (double ele 10 or single ele 20)
      hltTree->GetEntry(i);
      if( !HLT_SingleEle20) continue;

      //grab the rest of the important electron information 
      timer.StartSplit("Loading electron stuff");
      areAllEleBranchesOn = true;
      eTree->SetBranchStatus("*",1);
      eTree->GetEntry(i);     
 
      timer.StartSplit("Electron Cuts");
      //make a list of electrons passing our cuts
      std::vector< int > goodElectrons; 
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        //correct Pt
        if(!isMC) elePt->at(j) = energyScale.correctPt( elePt->at(j), eleSCEta->at(j), hiBin);
        if(isMC) elePt->at(j) = energyScaleMC.correctPt( elePt->at(j), eleSCEta->at(j), hiBin);
        
        //kinematic cuts
        if(elePt->at(j)< s.minElectronPt) continue;        
        if(TMath::Abs(eleSCEta->at(j)) > s.maxZRapEle) continue;
      
        //veto on transition region
        if(TMath::Abs(eleSCEta->at(j)) > 1.442 && TMath::Abs(eleSCEta->at(j)) < 1.556 ) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        //this dEta variable is the same as eledEtaSeedAtVtx (checked CMSSW)
        float dEta = TMath::Abs( eledEtaAtVtx->at(j) - eleSCEta->at(j) + eleSeedEta->at(j) );
        dEtaVsdEtaSeed->Fill(eledEtaAtVtx->at(j), dEta);
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverEBc->at(j), eleEoverPInv->at(j), eleIP3D->at(j))) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      bool moreThan2 = false;
      if(goodElectrons.size()>2) moreThan2 = true;

      //need to check if it is a tau daughter
      bool isTau = false;
      if(isMC){
        for(int j = 0; j<nMC; j++){
          //break out if you find a Z->tautau
          if( mcGMomPID->at(j) == 23 && TMath::Abs(mcMomPID->at(j)) == 15){
            isTau = true;
            break; 
          }
          //only look for electrons coming directly from a Z
          if( TMath::Abs(mcMomPID->at(j)) != 23) continue;  
          //if it decays to muons or electrons, break out
          if( TMath::Abs(mcPID->at(j)) == 13 || TMath::Abs(mcPID->at(j) == 11)) break;       
        }
      }

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

          //L1 trigger matching (1 L1 EG > 15 GeV)
          //bool isFirstElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 15.0);
          //bool isSecondElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 15.0);
          //if(! (isFirstElectronL1Matched || isSecondElectronL1Matched)) continue;

          //HLT trigger matching (1 HLT match > 20 GeV)
          bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
          bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
          if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
   
          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          double efficiencyArray[5] = {0};
          efficiencyArray[0] = zEff.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
          efficiencyArray[1] = zEffU.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
          efficiencyArray[2] = zEffD.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
          efficiencyArray[3] = zEffphotonU.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
          efficiencyArray[4] = zEffphotonD.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
          double efficiency = zEff.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
        
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " isOS? " << (int)isOppositeSign << std::endl;
          if( isOppositeSign){
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBinZDC,k)){
                if((isMC && !isTau) || !isMC){
                  massPeakOS[k]->Fill( Zcand.M() );
                  massPeakOS_withEff[k]->Fill( Zcand.M(), eventWeight/efficiency );
                  for( int l = 0; l< (isMC ? 1 : 5); l++){
                    pTOS_withEff[k][l]->Fill( Zcand.Pt(), 1.0/efficiencyArray[l] * eventWeight );
                    yOS_withEff[k][l]->Fill( Zcand.Rapidity(), 1.0/efficiencyArray[l] * eventWeight );
                    yieldOS_withEff[k][l]->Fill( 0.5, 1.0/efficiencyArray[l] * eventWeight );
                  }
  
                  lepPt->Fill( elePt->at(goodElectrons.at(j)) ,eventWeight);
                  lepPt->Fill( elePt->at(goodElectrons.at(j2)) ,eventWeight);
                  lepEta->Fill( eleEta->at(goodElectrons.at(j)) ,eventWeight);
                  lepEta->Fill( eleEta->at(goodElectrons.at(j2)) ,eventWeight);
                  lepPhi->Fill( elePhi->at(goodElectrons.at(j)) ,eventWeight);
                  lepPhi->Fill( elePhi->at(goodElectrons.at(j2)) ,eventWeight);
              
                  yields->Fill(k,eventWeight/efficiency);
                  massVsPt[k]->Fill(Zcand.M(), Zcand.Pt()); 
                  candPt[k]->Fill(Zcand.Pt());
                  candPt_unnormalized[k]->Fill(Zcand.Pt());
                  candPtFine[k]->Fill(Zcand.Pt());
                  candPtFiner[k]->Fill(Zcand.Pt());
                  candEta[k]->Fill(Zcand.Eta());
                  candY[k]->Fill(Zcand.Rapidity());
                  candPhi[k]->Fill(Zcand.Phi());
                  candPtVsM[k]->Fill(Zcand.Pt(),Zcand.M());
                  float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( elePhi->at(goodElectrons.at(j2)) - elePhi->at(goodElectrons.at(j)) )))/TMath::Pi(); 
                  candAcoVsM[k]->Fill(acoplanarity ,Zcand.M());
                  candAcoVsPt[k]->Fill(acoplanarity ,Zcand.Pt());
                  candAco[k]->Fill(acoplanarity);
                
                  //we have 3 different aco cuts for variations
                  bool passesAco[3] = {1 , 1, 1};
                  if(Zcand.Pt() < s.minPtCutForPhotons && acoplanarity < s.acoCutForPhotons) passesAco[0] = false;
                  if(Zcand.Pt() < s.minPtCutForPhotonsU && acoplanarity < s.acoCutForPhotonsU) passesAco[1] = false;
                  if(Zcand.Pt() < s.minPtCutForPhotonsD && acoplanarity < s.acoCutForPhotonsD) passesAco[2] = false;
                
                  if(Zcand.Pt() < s.minPtCutForPhotons){
                    massPeakOS_ptLT0p5_withEff[k]->Fill( Zcand.M(), 1.0/efficiency * eventWeight );
                  }
 
                  if( !passesAco[0]) massPeakOS_ptLT0p5acoLT0p001_withEff[k]->Fill( Zcand.M(), 1.0/efficiency*eventWeight);
                
                  for( int l = 0; l< (isMC ? 1 : 5); l++){
                    if( (!passesAco[0] && l<=2) || (!passesAco[1] && l==3) || (!passesAco[2] && l==4)){
                      pTOS_ptLT0p5acoLT0p001_withEff[k][l]->Fill( Zcand.Pt(), 1.0/efficiencyArray[l] * eventWeight );
                      yOS_ptLT0p5acoLT0p001_withEff[k][l]->Fill( Zcand.Rapidity(), 1.0/efficiencyArray[l]*eventWeight);
                      yieldOS_ptLT0p5acoLT0p001_withEff[k][l]->Fill( 0.5, 1.0/efficiencyArray[l]*eventWeight);
                    }
                  }     
                }//end MC if
                if( isMC && isTau ){
                  massPeakOS_TauTau_withEff[k]->Fill( Zcand.M(), 1.0/efficiency * eventWeight );
                  pTOS_TauTau_withEff[k]->Fill( Zcand.Pt(), 1.0/efficiency * eventWeight );
                  yOS_TauTau_withEff[k]->Fill( Zcand.Rapidity(), 1.0/efficiency * eventWeight );
                  yieldOS_TauTau_withEff[k]->Fill( 0.5, 1.0/efficiency * eventWeight );
                }
              }//end isInside if
            }//end for loop
          }else{
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBinZDC,k)){
                if((isMC && !isTau) || !isMC){
                  massPeakSS[k]->Fill( Zcand.M() );
                  massPeakSS_withEff[k]->Fill( Zcand.M(), eventWeight/efficiency );
                  for( int l = 0; l< (isMC ? 1 : 5); l++){
                    pTSS_withEff[k][l]->Fill( Zcand.Pt(), 1.0/efficiencyArray[l] * eventWeight );
                    ySS_withEff[k][l]->Fill( Zcand.Rapidity() , 1.0/efficiencyArray[l] * eventWeight);
                    yieldSS_withEff[k][l]->Fill( 0.5 , 1.0/efficiencyArray[l] * eventWeight);
                  }
                  yieldsSS->Fill(k,eventWeight/efficiency);
                
                  float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( elePhi->at(goodElectrons.at(j2)) - elePhi->at(goodElectrons.at(j)) )))/TMath::Pi(); 
                  if(Zcand.Pt() < s.minPtCutForPhotons){
                    massPeakSS_ptLT0p5_withEff[k]->Fill( Zcand.M(), 1.0/efficiency * eventWeight );
                    if( acoplanarity < s.acoCutForPhotons) massPeakSS_ptLT0p5acoLT0p001_withEff[k]->Fill( Zcand.M(), 1.0/efficiency*eventWeight);
                  }
                }//end MC if
                if( isMC && isTau ){
                  massPeakSS_TauTau_withEff[k]->Fill( Zcand.M(), 1.0/efficiency * eventWeight );
                  pTSS_TauTau_withEff[k]->Fill( Zcand.Pt(), 1.0/efficiency * eventWeight );
                  ySS_TauTau_withEff[k]->Fill( Zcand.Rapidity(), 1.0/efficiency * eventWeight );
                  yieldSS_TauTau_withEff[k]->Fill( 0.5, 1.0/efficiency * eventWeight );
                }     
              }
            }
          }
          if(isMC) continue;       
 
          //v2 calcuation
          if( isOppositeSign){
            //reference Q vectors
            TComplex Qp = TComplex(hiQVecMag[7], hiQVecAngle[7], true);
            TComplex Qn = TComplex(hiQVecMag[6], hiQVecAngle[6], true);
            TComplex Qmid = TComplex(hiQVecMag[9], hiQVecAngle[9], true);

            //signal Q vectors
            TComplex candQ = TComplex(1, 2*Zcand.Phi(), true);
            TComplex ele1Q = TComplex(1, 2*elePhi->at(goodElectrons.at(j)), true);
            TComplex ele2Q = TComplex(1, 2*elePhi->at(goodElectrons.at(j2)), true);

            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBinZDC,k)){
                //see equation 1 in HIN-16-007
                //'a' is Q1 and 'b' is Q2, c is Qmid
                TComplex Q1 = Qp;
                TComplex Q2 = Qn;
                if(Zcand.Eta()>0){
                  Q1 = Qn;
                  Q2 = Qp;
                }

                float num = (candQ*TComplex::Conjugate(Q1)).Re();
                float denom = (Q1*TComplex::Conjugate(Q2)).Re();
                float q1AndMid = (Q1*TComplex::Conjugate(Qmid)).Re();
                float q2AndMid = (Q2*TComplex::Conjugate(Qmid)).Re();
                v2Num[k]->Fill(0.5,num/efficiency);
                v2Denom[k]->Fill(0.5,denom);
                v2Q1Mid[k]->Fill(0.5,q1AndMid);
                v2Q2Mid[k]->Fill(0.5,q2AndMid);

                for(int l = 0; l<3; l++){
                  v2NumVsCent[l]->Fill(k,num/efficiencyArray[l]);
                  v2AvgEffVsCent[l]->Fill(k,1.0/efficiencyArray[l]);
                }
                v2DenomVsCent->Fill(k,denom);
                v2Q1MidVsCent->Fill(k,q1AndMid);
                v2Q2MidVsCent->Fill(k,q2AndMid);

                //electrons
                Q1 = Qp;
                Q2 = Qn;
                if(eleEta->at(goodElectrons.at(j))>0){
                  Q1 = Qn; 
                  Q2 = Qp;
                }
                float numEle1 = (ele1Q*TComplex::Conjugate(Q1)).Re();
                float denomEle1 = (Q1*TComplex::Conjugate(Q2)).Re();
                float q1AndMidEle1 = (Q1*TComplex::Conjugate(Qmid)).Re();
                float q2AndMidEle1 = (Q2*TComplex::Conjugate(Qmid)).Re();
                v2EleNum[k]->Fill(0.5,numEle1);
                v2EleDenom[k]->Fill(0.5,denomEle1);
                v2EleQ1Mid[k]->Fill(0.5,q1AndMidEle1);
                v2EleQ2Mid[k]->Fill(0.5,q2AndMidEle1);
              
                Q1 = Qp;
                Q2 = Qn;
                if(eleEta->at(goodElectrons.at(j2))>0){
                  Q1 = Qn; 
                  Q2 = Qp;
                }
                float numEle2 = (ele2Q*TComplex::Conjugate(Q1)).Re();
                float denomEle2 = (Q1*TComplex::Conjugate(Q2)).Re();
                float q1AndMidEle2 = (Q1*TComplex::Conjugate(Qmid)).Re();
                float q2AndMidEle2 = (Q2*TComplex::Conjugate(Qmid)).Re();
                v2EleNum[k]->Fill(0.5,numEle2);
                v2EleDenom[k]->Fill(0.5,denomEle2);
                v2EleQ1Mid[k]->Fill(0.5,q1AndMidEle2);
                v2EleQ2Mid[k]->Fill(0.5,q2AndMidEle2);

                v2EleNumVsCent->Fill(k,numEle1);
                v2EleNumVsCent->Fill(k,numEle2);
                v2EleDenomVsCent->Fill(k,denomEle1);
                v2EleDenomVsCent->Fill(k,denomEle2);
                v2EleQ1MidVsCent->Fill(k,q1AndMidEle1);
                v2EleQ1MidVsCent->Fill(k,q1AndMidEle2);
                v2EleQ2MidVsCent->Fill(k,q2AndMidEle1);
                v2EleQ2MidVsCent->Fill(k,q2AndMidEle2);
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
  
  for(int i = 0; i<nBins; i++){
    h.makeDifferential( candPt[i]);
    h.makeDifferential( candEta[i]);
    h.makeDifferential( candY[i]);
  
    for(int j = 0; j< (isMC ? 1 : 5); j++){
      h.makeDifferential( pTOS_withEff[i][j]);
      h.makeDifferential( pTSS_withEff[i][j]);
      h.makeDifferential( pTOS_ptLT0p5acoLT0p001_withEff[i][j]);
    }
    h.makeDifferential( pTOS_TauTau_withEff[i]);
    h.makeDifferential( pTSS_TauTau_withEff[i]);

    for(int j = 0; j< (isMC ? 1 : 5); j++){
      h.makeDifferential( yOS_withEff[i][j]);
      h.makeDifferential( ySS_withEff[i][j]);
      h.makeDifferential( yOS_ptLT0p5acoLT0p001_withEff[i][j]);
    }
    h.makeDifferential( yOS_TauTau_withEff[i]);
    h.makeDifferential( ySS_TauTau_withEff[i]);
  }
  

  TFile * output;
  if(!isMC){
    if(!doZDC) output = new TFile(Form("unmergedOutputs/Z2ee_%s_hiBin%d_%d_%d.root",outputTag.c_str(), hiBinVar, (int)isMC, jobNumber),"recreate");
    else       output = new TFile(Form("unmergedOutputs/Z2ee_%s_hiBinZDC%d_%d_%d.root",outputTag.c_str(), hiBinVar, (int)isMC, jobNumber),"recreate");
  }else{
    output = new TFile(Form("unmergedOutputs/Z2ee_%s_%d_%d.root",outputTag.c_str(), (int)isMC, jobNumber),"recreate");
  }
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    massPeakOS_withEff[i]->Write();
    massPeakSS_withEff[i]->Write();
    massPeakOS_ptLT0p5_withEff[i]->Write();
    massPeakSS_ptLT0p5_withEff[i]->Write();
    massPeakOS_ptLT0p5acoLT0p001_withEff[i]->Write();
    massPeakSS_ptLT0p5acoLT0p001_withEff[i]->Write();
    massPeakOS_TauTau_withEff[i]->Write();
    massPeakSS_TauTau_withEff[i]->Write();
    
    for(int j = 0; j< (isMC ? 1 : 5); j++){
      pTOS_withEff[i][j]->Write();
      pTSS_withEff[i][j]->Write();
      pTOS_ptLT0p5acoLT0p001_withEff[i][j]->Write(0);
    }
    pTOS_TauTau_withEff[i]->Write();
    pTSS_TauTau_withEff[i]->Write();
    
    for(int j = 0; j< (isMC ? 1 : 5); j++){
      yOS_withEff[i][j]->Write();
      ySS_withEff[i][j]->Write();
      yOS_ptLT0p5acoLT0p001_withEff[i][j]->Write(0);
    }
    yOS_TauTau_withEff[i]->Write();
    ySS_TauTau_withEff[i]->Write();
    
    for(int j = 0; j< (isMC ? 1 : 5); j++){
      yieldOS_withEff[i][j]->Write();
      yieldSS_withEff[i][j]->Write();
      yieldOS_ptLT0p5acoLT0p001_withEff[i][j]->Write(0);
    }
    yieldOS_TauTau_withEff[i]->Write();
    yieldSS_TauTau_withEff[i]->Write();

    massVsPt[i]->Write();
    candPt[i]->Write();
    candPt_unnormalized[i]->Write();
    candPtFine[i]->Write();
    candPtFiner[i]->Write();
    candEta[i]->Write();
    candY[i]->Write();
    candPhi[i]->Write();   
    candPtVsM[i]->Write();
    candAcoVsM[i]->Write();
    candAcoVsPt[i]->Write();
    candAco[i]->Write(); 

    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2Q1Mid[i]->Write();
    v2Q2Mid[i]->Write();
    v2EleNum[i]->Write();
    v2EleDenom[i]->Write();
    v2EleQ1Mid[i]->Write();
    v2EleQ2Mid[i]->Write();
  }
  nEvents->Write();
  yields->Write();
  yieldsSS->Write();
  for(int l = 0; l<3; l++){
    v2NumVsCent[l]->Write();
    v2AvgEffVsCent[l]->Write();
  }
  v2DenomVsCent->Write();
  v2Q1MidVsCent->Write();
  v2Q2MidVsCent->Write();
  v2EleNumVsCent->Write();
  v2EleDenomVsCent->Write();
  v2EleQ1MidVsCent->Write();
  v2EleQ2MidVsCent->Write();

  dEtaVsdEtaSeed->Write();
  
  lepPt->Write();
  lepEta->Write();
  lepPhi->Write();
    
  output->Close();

  timer.Stop();
  timer.Report();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 6)
  {
    std::cout << "Usage: Z_EE_Channel <fileList> <outputTag> <isMC> <job #> <total number of jobs>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string outputTag = argv[2];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());
  bool isMC = (bool)std::atoi(argv[3]);

  int job = std::atoi(argv[4]);
  int totalJobs = std::atoi(argv[5]);

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
   
  doZ2EE(listOfFiles, job, isMC, outputTag, 0);
  doZ2EE(listOfFiles, job, isMC, outputTag, 0, true);
  if(!isMC){
    doZ2EE(listOfFiles, job, isMC, outputTag, 1);
    doZ2EE(listOfFiles, job, isMC, outputTag, 2);
  }
  return 0; 
}

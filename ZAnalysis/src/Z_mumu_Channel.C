#include "include/HistNameHelper.h"
#include "include/MCReweight.h"
#include "include/centralityBin.h"
#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/ZEfficiency.h"
#include "include/HistNameHelper.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2mumu(std::vector< std::string > files, float etaCut, bool isMC, Settings s, std::string outFile, bool isTest, int job, int nJobs){
  TH1::SetDefaultSumw2();
  ZEfficiency zEff = ZEfficiency("resources/Z2mumu_Efficiencies.root", isMC);
  MCReweight * vzRW;
  if(isMC) vzRW = new MCReweight("resources/vzReweight.root");

  HistNameHelper h = HistNameHelper();
  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * nEvents = new TH1D("nEvents","nEvents",1,0,1);

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TH1D * massPeakOS_withEff[nBins]; 
  TH1D * massPeakSS_withEff[nBins]; 

  TH1D * massPeakOS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * massPeakSS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * massPeakOS_ptLT0p5_withEff[nBins]; 
  TH1D * massPeakSS_ptLT0p5_withEff[nBins]; 
  TH1D * massPeakOS_TauTau_withEff[nBins]; 
  
  TH1D * pTOS_withEff[nBins]; 
  TH1D * pTSS_withEff[nBins]; 
  TH1D * pTOS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * pTOS_TauTau_withEff[nBins]; 
  
  TH1D * yOS_withEff[nBins]; 
  TH1D * ySS_withEff[nBins]; 
  TH1D * yOS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * yOS_TauTau_withEff[nBins]; 
  
  TH1D * yieldOS_withEff[nBins]; 
  TH1D * yieldSS_withEff[nBins]; 
  TH1D * yieldOS_ptLT0p5acoLT0p001_withEff[nBins]; 
  TH1D * yieldOS_TauTau_withEff[nBins]; 

  TH1D * yields;
  TH1D * yieldsSS;
  TH1D * yields_TauTau;
  TProfile * v2Num[nBins];
  TProfile * v2NumVsCent;
  TProfile * v2Denom[nBins];
  TProfile * v2DenomVsCent;
  TProfile * v2Q1Mid[nBins];
  TProfile * v2Q1MidVsCent;
  TProfile * v2Q2Mid[nBins];
  TProfile * v2Q2MidVsCent;
  TProfile * v2AvgEff[nBins];
  TProfile * v2AvgEffVsCent;
  
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
  
  
  TProfile * v2MuNum[nBins];
  TProfile * v2MuNumVsCent;
  TProfile * v2MuDenom[nBins];
  TProfile * v2MuDenomVsCent;
  TProfile * v2MuQ1Mid[nBins];
  TProfile * v2MuQ1MidVsCent;
  TProfile * v2MuQ2Mid[nBins];
  TProfile * v2MuQ2MidVsCent;

  TProfile * avgMassVsPt[nBins];

  int hiBin;

  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_withEff[i] = new TH1D(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_withEff[i] = new TH1D(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    massPeakOS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("massPeakOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("massPeakSS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_ptLT0p5_withEff[i] = new TH1D(Form("massPeakOS_ptLT0p5_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_ptLT0p5_withEff[i] = new TH1D(Form("massPeakSS_ptLT0p5_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_TauTau_withEff[i] = new TH1D(Form("massPeakOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    pTOS_withEff[i] = new TH1D(Form("pTOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    pTSS_withEff[i] = new TH1D(Form("pTSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    pTOS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("pTOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    pTOS_TauTau_withEff[i] = new TH1D(Form("pTOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);

    int nRapBins = 16;
    if(etaCut < 2.11) nRapBins = 14;
    yOS_withEff[i] = new TH1D(Form("yOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-etaCut,etaCut);
    ySS_withEff[i] = new TH1D(Form("ySS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-etaCut,etaCut);
    yOS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("yOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-etaCut,etaCut);
    yOS_TauTau_withEff[i] = new TH1D(Form("yOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";y",nRapBins,-etaCut,etaCut);
    
    yieldOS_withEff[i] = new TH1D(Form("yieldOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
    yieldSS_withEff[i] = new TH1D(Form("yieldSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",1,0,1);
    yieldOS_ptLT0p5acoLT0p001_withEff[i] = new TH1D(Form("yieldOS_ptLT0p5acoLT0p001_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
    yieldOS_TauTau_withEff[i] = new TH1D(Form("yieldOS_TauTau_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",1,0,1);
    
    candPt[i] = new TH1D(Form("candPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candPtFine[i] = new TH1D(Form("candPtFine_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",400,0,200);
    candPtFiner[i] = new TH1D(Form("candPtFiner_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",800,0,200);
    candPt_unnormalized[i] = new TH1D(Form("candPt_unnormalized_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candEta[i] = new TH1D(Form("candEta_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-s.maxZRap,s.maxZRap);
    candY[i] = new TH1D(Form("candY_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",24,-s.maxZRap,s.maxZRap);
    candPhi[i] = new TH1D(Form("candPhi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-TMath::Pi(),TMath::Pi());
    candPtVsM[i] = new TH2D(Form("candPtVsM_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",40,0,10,60,60,120);
    candAcoVsM[i] = new TH2D(Form("candAcoVsM_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",50,0,0.01,60,60,120); 
    candAcoVsPt[i] = new TH2D(Form("candAcoVsPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",50,0,0.01,50,0,5); 
 
    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q1Mid[i] = new TProfile(Form("v2Q1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q2Mid[i] = new TProfile(Form("v2Q2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2AvgEff[i] = new TProfile(Form("v2AvgEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
 
    v2MuNum[i] = new TProfile(Form("v2MuNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuDenom[i] = new TProfile(Form("v2MuDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuQ1Mid[i] = new TProfile(Form("v2MuQ1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuQ2Mid[i] = new TProfile(Form("v2MuQ2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
   
    avgMassVsPt[i] = new TProfile(Form("avgMassVsPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",10,0,5);
  }  
  v2NumVsCent = new TProfile("v2NumVsCent","",nBins,0,nBins);
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2Q1MidVsCent = new TProfile("v2Q1MidVsCent","",nBins,0,nBins); 
  v2Q2MidVsCent = new TProfile("v2Q2MidVsCent","",nBins,0,nBins); 
  v2AvgEffVsCent = new TProfile("v2AvgEffVsCent","",nBins,0,nBins);
  yields = new TH1D("yields","yields",nBins,0,nBins);   
  yieldsSS = new TH1D("yieldsSS","yieldsSS",nBins,0,nBins);   
  yields_TauTau = new TH1D("yields_TauTau","yields_TauTau",nBins,0,nBins);   
 
  v2MuNumVsCent = new TProfile("v2MuNumVsCent","",nBins,0,nBins);
  v2MuDenomVsCent = new TProfile("v2MuDenomVsCent","",nBins,0,nBins);
  v2MuQ1MidVsCent = new TProfile("v2MuQ1MidVsCent","",nBins,0,nBins); 
  v2MuQ2MidVsCent = new TProfile("v2MuQ2MidVsCent","",nBins,0,nBins); 

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    std::string tree = ""; 
    if(isMC) tree = "dimucontana_mc";
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),tree);
    unsigned int maxEvt = isTest ? 100000 : v.GetEntries(); 
    for(unsigned int i = 0; i<maxEvt; i++){
      if( i%nJobs != (unsigned int)job) continue;
      v.GetEntry(i);
      if(i%5000==0) std::cout << i << "/" << maxEvt << std::endl;

      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      hiBin = cb.getHiBinFromhiHFSides(v.HFsumETPlus(), v.HFsumETMinus(), (isMC ? 3 : 0));
      float eventWeight = 1.0;
      if(isMC) eventWeight = vzRW->reweightFactor( v.bestvtxZ() ) * c.findNcoll( hiBin );
      nEvents->Fill(0.5,eventWeight);
      
      //check out trigger 
      if( !(v.trigHLT()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12]) ) continue;

      //tag DY products in our MC sample
      bool isTau = false;
      if(isMC){
        for(unsigned int k = 0; k<v.candSize_gen(); k++){
          if( TMath::Abs(v.DecayID_gen()[k]) == 23 ){
            isTau = false;
            break;
          }
          if( TMath::Abs(v.DecayID_gen()[k]) == 15 ){
            isTau = true;
            break;
          }
        }
      }

      for(unsigned int j = 0; j<v.candSize(); j++){
        if( v.mass()[j] < s.zMassRange[0] || v.mass()[j] > s.zMassRange[1]) continue; 
        if( !(v.pTD1()[j] > s.minMuonPt )) continue;
        if( !(v.pTD2()[j] > s.minMuonPt )) continue;
        if( TMath::Abs(v.EtaD1()[j]) > etaCut ) continue;
        if( TMath::Abs(v.EtaD2()[j]) > etaCut ) continue;
        if( !(v.tightCand(j,"POG"))) continue;//tight Muon 1 && tight Muon 2      
        if( !(v.VtxProb()[j] >0.001)) continue; 
 
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[j] != v.chargeD2()[j];
        double efficiency = zEff.getEfficiency( v.y()[j], v.pT()[j] , hiBin );

        if( isOppositeSign){
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin( hiBin ,k)){
    
              if((isMC && !isTau) || !isMC){
                massPeakOS[k]->Fill( v.mass()[j] );
                massPeakOS_withEff[k]->Fill( v.mass()[j], 1.0/efficiency * eventWeight );
                pTOS_withEff[k]->Fill( v.pT()[j], 1.0/efficiency * eventWeight );
                yOS_withEff[k]->Fill( v.y()[j], 1.0/efficiency * eventWeight );
                yieldOS_withEff[k]->Fill( 0.5, 1.0/efficiency * eventWeight );
                
                yields->Fill(k,1.0/efficiency * eventWeight);
                candPt[k]->Fill(v.pT()[j]);
                candPtFine[k]->Fill(v.pT()[j]);
                candPtFiner[k]->Fill(v.pT()[j]);
                candPt_unnormalized[k]->Fill(v.pT()[j]);
                candEta[k]->Fill(v.eta()[j]);
                candY[k]->Fill(v.y()[j]);
                candPhi[k]->Fill(v.phi()[j]);
                avgMassVsPt[k]->Fill(v.pT()[j], v.mass()[j]);
                candPtVsM[k]->Fill(v.pT()[j],v.mass()[j]);
                float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( v.PhiD1()[j] - v.PhiD2()[j] )))/TMath::Pi(); 
                candAcoVsM[k]->Fill(acoplanarity ,v.mass()[j]);
                candAcoVsPt[k]->Fill(acoplanarity ,v.pT()[j]);
                
                if(v.pT()[j] < s.minPtCutForPhotons){
                  massPeakOS_ptLT0p5_withEff[k]->Fill( v.mass()[j], 1.0/efficiency * eventWeight );
                  if( acoplanarity < s.acoCutForPhotons){
                    massPeakOS_ptLT0p5acoLT0p001_withEff[k]->Fill( v.mass()[j], 1.0/efficiency*eventWeight);
                    pTOS_ptLT0p5acoLT0p001_withEff[k]->Fill( v.pT()[j], 1.0/efficiency*eventWeight);
                    yOS_ptLT0p5acoLT0p001_withEff[k]->Fill( v.y()[j], 1.0/efficiency*eventWeight);
                    yieldOS_ptLT0p5acoLT0p001_withEff[k]->Fill( 0.5, 1.0/efficiency*eventWeight);
                  }
                }
              }
              if( isMC && isTau ){
                massPeakOS_TauTau_withEff[k]->Fill( v.mass()[j], 1.0/efficiency * eventWeight );
                pTOS_TauTau_withEff[k]->Fill( v.pT()[j], 1.0/efficiency * eventWeight );
                yOS_TauTau_withEff[k]->Fill( v.y()[j], 1.0/efficiency * eventWeight );
                yieldOS_TauTau_withEff[k]->Fill( 0.5, 1.0/efficiency * eventWeight );
                yields_TauTau->Fill(k,1.0/efficiency * eventWeight);
              }
            }
          }
        }else{
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin( hiBin ,k)){
              if((isMC && !isTau) || !isMC){
                candPtFine[k]->Fill(v.pT()[j],-1);
                candPtFiner[k]->Fill(v.pT()[j],-1);
                massPeakSS[k]->Fill( v.mass()[j] );
                massPeakSS_withEff[k]->Fill( v.mass()[j] , 1.0/efficiency * eventWeight);
                pTSS_withEff[k]->Fill( v.pT()[j] , 1.0/efficiency * eventWeight);
                ySS_withEff[k]->Fill( v.y()[j] , 1.0/efficiency * eventWeight);
                yieldSS_withEff[k]->Fill( 0.5 , 1.0/efficiency * eventWeight);
                yieldsSS->Fill(k,1.0/efficiency * eventWeight);
                
                float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( v.PhiD1()[j] - v.PhiD2()[j] )))/TMath::Pi(); 
                if(v.pT()[j] < s.minPtCutForPhotons){
                  massPeakSS_ptLT0p5_withEff[k]->Fill( v.mass()[j], 1.0/efficiency * eventWeight );
                  if( acoplanarity < s.acoCutForPhotons) massPeakSS_ptLT0p5acoLT0p001_withEff[k]->Fill( v.mass()[j], 1.0/efficiency*eventWeight);
                }
              }
            }
          }
        }
        if(isMC) continue;

        //v2 calcuation
        TComplex Qp = TComplex(v.ephfpQ()[1], v.ephfpAngle()[1], true);
        TComplex Qn = TComplex(v.ephfmQ()[1], v.ephfmAngle()[1], true);
        TComplex Qmid = TComplex(v.eptrackmidQ()[1],v.eptrackmidAngle()[1],true);

        TComplex candQ = TComplex(1, 2*v.phi()[j], true);
        TComplex mu1Q = TComplex(1, 2*v.PhiD1()[j], true);
        TComplex mu2Q = TComplex(1, 2*v.PhiD2()[j], true);
        if( isOppositeSign){
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin( hiBin ,k)){
              TComplex Q1 = Qp;
              TComplex Q2 = Qn;
              if(v.eta()[j]>0){
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
              v2AvgEff[k]->Fill(0.5,1.0/efficiency);

              v2NumVsCent->Fill(k,num/efficiency);
              v2DenomVsCent->Fill(k,denom);
              v2Q1MidVsCent->Fill(k,q1AndMid);
              v2Q2MidVsCent->Fill(k,q2AndMid);
              v2AvgEffVsCent->Fill(k,1.0/efficiency);

              //muons
              Q1 = Qp;
              Q2 = Qn;
              if(v.EtaD1()[j]>0){
                Q1 = Qn; 
                Q2 = Qp;
              }
              float numMu1 = (mu1Q*TComplex::Conjugate(Q1)).Re();
              float denomMu1 = (Q1*TComplex::Conjugate(Q2)).Re();
              float q1AndMidMu1 = (Q1*TComplex::Conjugate(Qmid)).Re();
              float q2AndMidMu1 = (Q2*TComplex::Conjugate(Qmid)).Re();
              v2MuNum[k]->Fill(0.5,numMu1);
              v2MuDenom[k]->Fill(0.5,denomMu1);
              v2MuQ1Mid[k]->Fill(0.5,q1AndMidMu1);
              v2MuQ2Mid[k]->Fill(0.5,q2AndMidMu1);
              
              Q1 = Qp;
              Q2 = Qn;
              if(v.EtaD2()[j]>0){
                Q1 = Qn; 
                Q2 = Qp;
              }
              float numMu2 = (mu2Q*TComplex::Conjugate(Q1)).Re();
              float denomMu2 = (Q1*TComplex::Conjugate(Q2)).Re();
              float q1AndMidMu2 = (Q1*TComplex::Conjugate(Qmid)).Re();
              float q2AndMidMu2 = (Q2*TComplex::Conjugate(Qmid)).Re();
              v2MuNum[k]->Fill(0.5,numMu2);
              v2MuDenom[k]->Fill(0.5,denomMu2);
              v2MuQ1Mid[k]->Fill(0.5,q1AndMidMu2);
              v2MuQ2Mid[k]->Fill(0.5,q2AndMidMu2);

              v2MuNumVsCent->Fill(k,numMu1);
              v2MuNumVsCent->Fill(k,numMu2);
              v2MuDenomVsCent->Fill(k,denomMu1);
              v2MuDenomVsCent->Fill(k,denomMu2);
              v2MuQ1MidVsCent->Fill(k,q1AndMidMu1);
              v2MuQ1MidVsCent->Fill(k,q1AndMidMu2);
              v2MuQ2MidVsCent->Fill(k,q2AndMidMu1);
              v2MuQ2MidVsCent->Fill(k,q2AndMidMu2);
            }
          }
        }
      }
    }
  }

  for(int i = 0; i<nBins; i++){
    h.makeDifferential( candPt[i]);
    h.makeDifferential( candEta[i]);
    h.makeDifferential( candY[i]);

    h.makeDifferential( pTOS_withEff[i]);
    h.makeDifferential( pTSS_withEff[i]);
    h.makeDifferential( pTOS_ptLT0p5acoLT0p001_withEff[i]);
    h.makeDifferential( pTOS_TauTau_withEff[i]);

    h.makeDifferential( yOS_withEff[i]);
    h.makeDifferential( ySS_withEff[i]);
    h.makeDifferential( yOS_ptLT0p5acoLT0p001_withEff[i]);
    h.makeDifferential( yOS_TauTau_withEff[i]);
  }

  TFile * output;
  if(!isTest){
    if(!isMC) output = new TFile(Form("Z2mumu_%d_%s_job%d.root",(int)(etaCut*10), outFile.c_str(), job),"recreate");
    if(isMC) output = new TFile(Form("Z2mumu_MC_%d_%s_job%d.root",(int)(etaCut*10), outFile.c_str(), job),"recreate");
  }else{
    if(!isMC) output = new TFile(Form("Z2mumu_%d_%s_job%d_TEST.root",(int)(etaCut*10), outFile.c_str(), job),"recreate");
    if(isMC) output = new TFile(Form("Z2mumu_MC_%d_%s_job%d_TEST.root",(int)(etaCut*10), outFile.c_str(), job),"recreate");
  }
  for(int i = 0; i<nBins; i++){
    float entries = candPtFine[i]->Integral();
    candPtFine[i]->Scale(1.0/entries);
    entries = candPtFiner[i]->Integral();
    candPtFiner[i]->Scale(1.0/entries);

    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    massPeakOS_withEff[i]->Write();
    massPeakSS_withEff[i]->Write();
    massPeakOS_ptLT0p5_withEff[i]->Write(0);
    massPeakSS_ptLT0p5_withEff[i]->Write(0);
    massPeakOS_ptLT0p5acoLT0p001_withEff[i]->Write(0);
    massPeakSS_ptLT0p5acoLT0p001_withEff[i]->Write(0);
    massPeakOS_TauTau_withEff[i]->Write();

    pTOS_withEff[i]->Write();
    pTSS_withEff[i]->Write();
    pTOS_ptLT0p5acoLT0p001_withEff[i]->Write(0);
    pTOS_TauTau_withEff[i]->Write();
    
    yOS_withEff[i]->Write();
    ySS_withEff[i]->Write();
    yOS_ptLT0p5acoLT0p001_withEff[i]->Write(0);
    yOS_TauTau_withEff[i]->Write();
    
    yieldOS_withEff[i]->Write();
    yieldSS_withEff[i]->Write();
    yieldOS_ptLT0p5acoLT0p001_withEff[i]->Write(0);
    yieldOS_TauTau_withEff[i]->Write();

    candPt[i]->Write();
    candPtFine[i]->Write();
    candPtFiner[i]->Write();
    candPt_unnormalized[i]->Write();
    candEta[i]->Write();
    candY[i]->Write();
    candPhi[i]->Write();  
    candPtVsM[i]->Write(); 
    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2Q1Mid[i]->Write();
    v2Q2Mid[i]->Write();
    v2AvgEff[i]->Write();
    v2MuNum[i]->Write();
    v2MuDenom[i]->Write();
    v2MuQ1Mid[i]->Write();
    v2MuQ2Mid[i]->Write();
    candAcoVsM[i]->Write();
    candAcoVsPt[i]->Write();
    avgMassVsPt[i]->Write();
  }
  
  v2NumVsCent->Write();
  v2DenomVsCent->Write();
  v2Q1MidVsCent->Write();
  v2Q2MidVsCent->Write();
  v2AvgEffVsCent->Write();
  v2MuNumVsCent->Write();
  v2MuDenomVsCent->Write();
  v2MuQ1MidVsCent->Write();
  v2MuQ2MidVsCent->Write();

  yields->Write();
  yieldsSS->Write();
  yields_TauTau->Write();
  nEvents->Write();

  output->Close();

  if(isMC) delete vzRW;

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 7)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList> <outFileTag> <isMC> < isTest> <jobNumber> <nJobs>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string outFile = argv[2];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());
 
  bool isMC = (bool)std::atoi(argv[3]);
  bool isTest = (bool)std::atoi(argv[4]);
  bool job = (bool)std::atoi(argv[5]);
  bool nJobs = (bool)std::atoi(argv[6]);

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
      listOfFiles.push_back(buffer);
      line++;
    }
  }

  Settings s = Settings();
   
  doZ2mumu(listOfFiles, s.maxZRap, isMC, s, outFile, isTest, job, nJobs);
  doZ2mumu(listOfFiles, s.maxZRapEle,isMC, s, outFile, isTest, job, nJobs);
  return 0; 
}

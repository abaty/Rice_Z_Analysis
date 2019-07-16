#include "include/centralityBin.h"
#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/ZEfficiency.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2mumu(std::vector< std::string > files, float etaCut, bool isMC, Settings s){
  TH1::SetDefaultSumw2();
  ZEfficiency zEff = ZEfficiency("resources/Z2mumu_Efficiencies.root");

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TH1D * massPeakOS_withEff[nBins]; 
  TH1D * massPeakSS_withEff[nBins]; 
  TH1D * yields;
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
  TH1D * candPt_unnormalized[nBins];
  TH1D * candEta[nBins];
  TH1D * candY[nBins];
  TH1D * candPhi[nBins];
  
  
  TProfile * v2MuNum[nBins];
  TProfile * v2MuNumVsCent;
  TProfile * v2MuDenom[nBins];
  TProfile * v2MuDenomVsCent;
  TProfile * v2MuQ1Mid[nBins];
  TProfile * v2MuQ1MidVsCent;
  TProfile * v2MuQ2Mid[nBins];
  TProfile * v2MuQ2MidVsCent;

  int hiBin;

  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_withEff[i] = new TH1D(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_withEff[i] = new TH1D(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    candPt[i] = new TH1D(Form("candPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candPt_unnormalized[i] = new TH1D(Form("candPt_unnormalized_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candEta[i] = new TH1D(Form("candEta_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-s.maxZRap,s.maxZRap);
    candY[i] = new TH1D(Form("candY_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",24,-s.maxZRap,s.maxZRap);
    candPhi[i] = new TH1D(Form("candPhi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-TMath::Pi(),TMath::Pi());
 
    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q1Mid[i] = new TProfile(Form("v2Q1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q2Mid[i] = new TProfile(Form("v2Q2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2AvgEff[i] = new TProfile(Form("v2AvgEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
 
    v2MuNum[i] = new TProfile(Form("v2MuNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuDenom[i] = new TProfile(Form("v2MuDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuQ1Mid[i] = new TProfile(Form("v2MuQ1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuQ2Mid[i] = new TProfile(Form("v2MuQ2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
  }  
  v2NumVsCent = new TProfile("v2NumVsCent","",nBins,0,nBins);
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2Q1MidVsCent = new TProfile("v2Q1MidVsCent","",nBins,0,nBins); 
  v2Q2MidVsCent = new TProfile("v2Q2MidVsCent","",nBins,0,nBins); 
  v2AvgEffVsCent = new TProfile("v2AvgEffVsCent","",nBins,0,nBins);
  yields = new TH1D("yields","yields",nBins,0,nBins);   
 
  v2MuNumVsCent = new TProfile("v2MuNumVsCent","",nBins,0,nBins);
  v2MuDenomVsCent = new TProfile("v2MuDenomVsCent","",nBins,0,nBins);
  v2MuQ1MidVsCent = new TProfile("v2MuQ1MidVsCent","",nBins,0,nBins); 
  v2MuQ2MidVsCent = new TProfile("v2MuQ2MidVsCent","",nBins,0,nBins); 

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),""); 
    for(unsigned int i = 0; i<v.GetEntries(); i++){
      v.GetEntry(i);

      //check out trigger 
      if( !(v.trigHLT()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12]) ) continue;

      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      hiBin = cb.getHiBinFromhiHFSides(v.HFsumETPlus(), v.HFsumETMinus(),0);

      if(i%1000==0) std::cout << i << std::endl;

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
              massPeakOS[k]->Fill( v.mass()[j] );
              massPeakOS_withEff[k]->Fill( v.mass()[j], 1.0/efficiency );
              yields->Fill(k,1.0/efficiency);
              candPt[k]->Fill(v.pT()[j]);
              candPt_unnormalized[k]->Fill(v.pT()[j]);
              candEta[k]->Fill(v.eta()[j]);
              candY[k]->Fill(v.y()[j]);
              candPhi[k]->Fill(v.phi()[j]);
            }
          }
        }else{
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin( hiBin ,k)){
              massPeakSS[k]->Fill( v.mass()[j] );
              massPeakSS_withEff[k]->Fill( v.mass()[j] , 1.0/efficiency);
              yields->Fill(k,-1.0/efficiency);
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
    //make differential quantities
    for(int j = 1; j<candPt[i]->GetNbinsX()+1; j++){
      candPt[i]->SetBinContent(j,candPt[i]->GetBinContent(j)/candPt[i]->GetBinWidth(j));
      candPt[i]->SetBinError(j,candPt[i]->GetBinError(j)/candPt[i]->GetBinWidth(j));
    }
    for(int j = 1; j<candEta[i]->GetNbinsX()+1; j++){
      candEta[i]->SetBinContent(j,candEta[i]->GetBinContent(j)/candEta[i]->GetBinWidth(j));
      candEta[i]->SetBinError(j,candEta[i]->GetBinError(j)/candEta[i]->GetBinWidth(j));
    }
    for(int j = 1; j<candY[i]->GetNbinsX()+1; j++){
      candY[i]->SetBinContent(j,candY[i]->GetBinContent(j)/candY[i]->GetBinWidth(j));
      candY[i]->SetBinError(j,candY[i]->GetBinError(j)/candY[i]->GetBinWidth(j));
    }

    massPeakOS[i]->SetDirectory(0);
    massPeakSS[i]->SetDirectory(0);
    massPeakOS_withEff[i]->SetDirectory(0);
    massPeakSS_withEff[i]->SetDirectory(0);
    candPt[i]->SetDirectory(0);
    candPt_unnormalized[i]->SetDirectory(0);
    candEta[i]->SetDirectory(0);
    candY[i]->SetDirectory(0);
    candPhi[i]->SetDirectory(0);
    v2Num[i]->SetDirectory(0);
    v2Denom[i]->SetDirectory(0);
    v2Q1Mid[i]->SetDirectory(0);
    v2Q2Mid[i]->SetDirectory(0);
    v2AvgEff[i]->SetDirectory(0);
    v2MuNum[i]->SetDirectory(0);
    v2MuDenom[i]->SetDirectory(0);
    v2MuQ1Mid[i]->SetDirectory(0);
    v2MuQ2Mid[i]->SetDirectory(0);
    yields->SetDirectory(0);
  }

  v2DenomVsCent->SetDirectory(0);
  v2NumVsCent->SetDirectory(0);
  v2Q1MidVsCent->SetDirectory(0);
  v2Q2MidVsCent->SetDirectory(0);
  v2AvgEffVsCent->SetDirectory(0); 
 
  v2MuDenomVsCent->SetDirectory(0);
  v2MuNumVsCent->SetDirectory(0);
  v2MuQ1MidVsCent->SetDirectory(0);
  v2MuQ2MidVsCent->SetDirectory(0);

  TFile * output;
  if(!isMC) output = new TFile(Form("Z2mumu_%d.root",(int)(etaCut*10)),"recreate");
  if(isMC) output = new TFile(Form("Z2mumu_MC_%d.root",(int)(etaCut*10)),"recreate");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    massPeakOS_withEff[i]->Write();
    massPeakSS_withEff[i]->Write();
    candPt[i]->Write();
    candPt_unnormalized[i]->Write();
    candEta[i]->Write();
    candY[i]->Write();
    candPhi[i]->Write();   
    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2Q1Mid[i]->Write();
    v2Q2Mid[i]->Write();
    v2AvgEff[i]->Write();
    v2MuNum[i]->Write();
    v2MuDenom[i]->Write();
    v2MuQ1Mid[i]->Write();
    v2MuQ2Mid[i]->Write();
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

  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList> <isMC>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());
 
  bool isMC = (bool)std::atoi(argv[2]);

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
   
  doZ2mumu(listOfFiles, s.maxZRap, isMC, s);
  doZ2mumu(listOfFiles, s.maxZRapEle,isMC, s);
  return 0; 
}

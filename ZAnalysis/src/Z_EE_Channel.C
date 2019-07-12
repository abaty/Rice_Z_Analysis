#include "include/electronEnergyScale.h"
#include "include/electronSelector.h"
#include "include/electronTriggerMatching.h"
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

void doZ2EE(std::vector< std::string > files, int jobNumber){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  //switch between single electron 20 and double electron 10
  bool doSingleEle20 = true;

  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronEnergyScale energyScale = ElectronEnergyScale("data");
  ZEfficiency zEff = ZEfficiency("resources/Z2ee_EfficiencyMC_0.root");
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TH1D * massPeakOS_withEff[nBins]; 
  TH1D * massPeakSS_withEff[nBins]; 
  TH2D * massVsPt[nBins];
  TH1D * candPt[nBins];
  TH1D * candEta[nBins];
  TH1D * candY[nBins];
  TH1D * candPhi[nBins];
  TH1D * yields;
  
  TProfile * v2Num[nBins];
  TProfile * v2NumVsCent;
  TProfile * v2Denom[nBins];
  TProfile * v2DenomVsCent;
  TProfile * v2Q1Mid[nBins];
  TProfile * v2Q1MidVsCent;
  TProfile * v2Q2Mid[nBins];
  TProfile * v2Q2MidVsCent;
  
  TProfile * v2EleNum[nBins];
  TProfile * v2EleNumVsCent;
  TProfile * v2EleDenom[nBins];
  TProfile * v2EleDenomVsCent;
  TProfile * v2EleQ1Mid[nBins];
  TProfile * v2EleQ1MidVsCent;
  TProfile * v2EleQ2Mid[nBins];
  TProfile * v2EleQ2MidVsCent;

  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakOS_withEff[i] = new TH1D(Form("massPeakOS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS_withEff[i] = new TH1D(Form("massPeakSS_withEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massVsPt[i] = new TH2D(Form("massVsPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}};p_{T}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1],s.nZPtBins-1,s.zPtBins);
    candPt[i] = new TH1D(Form("candPt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",s.nZPtBins-1,s.zPtBins);
    candEta[i] = new TH1D(Form("candEta_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-s.maxZRapEle,s.maxZRapEle);
    candY[i] = new TH1D(Form("candY_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-s.maxZRapEle,s.maxZRapEle);
    candPhi[i] = new TH1D(Form("candPhi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";p_{T}",20,-TMath::Pi(),TMath::Pi());

    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q1Mid[i] = new TProfile(Form("v2Q1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Q2Mid[i] = new TProfile(Form("v2Q2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
    v2EleNum[i] = new TProfile(Form("v2EleNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleDenom[i] = new TProfile(Form("v2EleDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleQ1Mid[i] = new TProfile(Form("v2EleQ1Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleQ2Mid[i] = new TProfile(Form("v2EleQ2Mid_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
  }  
  yields = new TH1D("yields","yields",nBins,0,nBins);  
  v2NumVsCent = new TProfile("v2NumVsCent","",nBins,0,nBins);
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2Q1MidVsCent = new TProfile("v2Q1MidVsCent","",nBins,0,nBins);
  v2Q2MidVsCent = new TProfile("v2Q2MidVsCent","",nBins,0,nBins);

  v2EleNumVsCent = new TProfile("v2EleNumVsCent","",nBins,0,nBins);
  v2EleDenomVsCent = new TProfile("v2EleDenomVsCent","",nBins,0,nBins);
  v2EleQ1MidVsCent = new TProfile("v2EleQ1MidVsCent","",nBins,0,nBins);
  v2EleQ2MidVsCent = new TProfile("v2EleQ2MidVsCent","",nBins,0,nBins);

  int nEle;
  int hiBin;
  float vz;

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
  std::vector< float > * eleSCPhi = 0;
  std::vector< float > * eleHoverEBc = 0;
  //std::vector< float > * eleD0 = 0;
  //std::vector< float > * eleDz = 0;
  std::vector< float > * eleIP3D = 0;
  std::vector< float > * eleEoverPInv = 0;

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
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverEBc",&eleHoverEBc);
    //eTree->SetBranchAddress("eleD0",&eleD0);
    //eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleIP3D",&eleIP3D);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("vz",&vz);
    evtTree->SetBranchAddress("hiNevtPlane",&hiNevtPlane);  
    evtTree->SetBranchAddress("hiQVecMag",hiQVecMag);  
    evtTree->SetBranchAddress("hiQVecAngle",hiQVecAngle);  
  
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

    TTree * L1Tree = (TTree*)in->Get("l1object/L1UpgradeFlatTree");
    L1Tree->SetBranchAddress("nEGs",&(eTrig.L1nEGs));
    L1Tree->SetBranchAddress("egEta", &(eTrig.L1egEta));
    L1Tree->SetBranchAddress("egPhi", &(eTrig.L1egPhi));
    L1Tree->SetBranchAddress("egEt", &(eTrig.L1egEt));

    TTree * HLTObjTree;
    if(!doSingleEle20){
      HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIDoubleEle10GsfMass50_v");
    } else {
      HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIEle20Gsf_v");
    }
    HLTObjTree->SetBranchAddress("eta",&(eTrig.HLTEta));
    HLTObjTree->SetBranchAddress("phi",&(eTrig.HLTPhi));
    HLTObjTree->SetBranchAddress("pt",&(eTrig.HLTPt));

    bool areAllEleBranchesOn = true;
    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
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
      if((!doSingleEle20) && (!HLT_DoubleEle10)) continue;
      if( doSingleEle20 && (!HLT_SingleEle20)) continue;

      timer.StartSplit("Checking Evt Selections");
      //event selections
      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;
      skimTree->GetEntry(i);
      if(! (pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter)) continue;

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
        elePt->at(j) = energyScale.correctPt( elePt->at(j), eleSCEta->at(j), hiBin);
        
        //kinematic cuts
        if(elePt->at(j)< s.minElectronPt) continue;        
        if(TMath::Abs(eleSCEta->at(j)) > s.maxZRapEle) continue;
      
        //veto on transition region
        if(TMath::Abs(eleSCEta->at(j)) > 1.442 && TMath::Abs(eleSCEta->at(j)) < 1.556 ) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        float dEta = TMath::Abs( eledEtaAtVtx->at(j) );
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverEBc->at(j), eleEoverPInv->at(j), eleIP3D->at(j))) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      bool moreThan2 = false;
      if(goodElectrons.size()>2) moreThan2 = true;

      timer.StartSplit("Loading HLT/L1 Object stuff");
      //get trigger matching stuff
      L1Tree->GetEntry(i);
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
          bool isFirstElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 15.0);
          bool isSecondElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 15.0);
          if(! (isFirstElectronL1Matched || isSecondElectronL1Matched)) continue;

          //HLT trigger matching (1 HLT match > 20 GeV)
          if(doSingleEle20){
            bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
            bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
            if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
          } else {  //(2 HLT matches > 10 GeV)
            bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 10.0);
            bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 10.0);
            if(! (isFirstElectronHLTMatched && isSecondElectronHLTMatched)) continue;
          }
   
          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          double efficiency = zEff.getEfficiency( Zcand.Rapidity(), Zcand.Pt() , hiBin );
        
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " isOS? " << (int)isOppositeSign << std::endl;
          if( isOppositeSign){
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBin,k)){
                massPeakOS[k]->Fill( Zcand.M() );
                massPeakOS_withEff[k]->Fill( Zcand.M(), 1.0/efficiency );
                yields->Fill(k,1.0/efficiency);
                massVsPt[k]->Fill(Zcand.M(), Zcand.Pt()); 
                candPt[k]->Fill(Zcand.Pt());
                candEta[k]->Fill(Zcand.Eta());
                candY[k]->Fill(Zcand.Rapidity());
                candPhi[k]->Fill(Zcand.Phi());
              }
            }
          }else{
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBin,k)){
                massPeakSS[k]->Fill( Zcand.M() );
                massPeakSS_withEff[k]->Fill( Zcand.M(), 1.0/efficiency );
                yields->Fill(k,-1.0/efficiency);
              }
            }
          }

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
              if(c.isInsideBin(hiBin,k)){
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

                v2NumVsCent->Fill(k,num);
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

  timer.StartSplit("End of analysis");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->SetDirectory(0);
    massPeakSS[i]->SetDirectory(0);
    massPeakOS_withEff[i]->SetDirectory(0);
    massPeakSS_withEff[i]->SetDirectory(0);
    massVsPt[i]->SetDirectory(0);
    candPt[i]->SetDirectory(0);
    candEta[i]->SetDirectory(0);
    candY[i]->SetDirectory(0);
    candPhi[i]->SetDirectory(0);

    v2Num[i]->SetDirectory(0);
    v2Denom[i]->SetDirectory(0);
    v2Q1Mid[i]->SetDirectory(0);
    v2Q2Mid[i]->SetDirectory(0);
    v2EleNum[i]->SetDirectory(0);
    v2EleDenom[i]->SetDirectory(0);
    v2EleQ1Mid[i]->SetDirectory(0);
    v2EleQ2Mid[i]->SetDirectory(0);
  }

  yields->SetDirectory(0);
  v2NumVsCent->SetDirectory(0);
  v2DenomVsCent->SetDirectory(0);
  v2Q1MidVsCent->SetDirectory(0);
  v2Q2MidVsCent->SetDirectory(0);

  v2EleNumVsCent->SetDirectory(0);
  v2EleDenomVsCent->SetDirectory(0);
  v2EleQ1MidVsCent->SetDirectory(0);
  v2EleQ2MidVsCent->SetDirectory(0);

  TFile * output = new TFile(Form("unmergedOutputs/Z2ee_%d.root",jobNumber),"recreate");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    massPeakOS_withEff[i]->Write();
    massPeakSS_withEff[i]->Write();
    massVsPt[i]->Write();
    candPt[i]->Write();
    candEta[i]->Write();
    candY[i]->Write();
    candPhi[i]->Write();   
 
    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2Q1Mid[i]->Write();
    v2Q2Mid[i]->Write();
    v2EleNum[i]->Write();
    v2EleDenom[i]->Write();
    v2EleQ1Mid[i]->Write();
    v2EleQ2Mid[i]->Write();
  }
  yields->Write();
  v2NumVsCent->Write();
  v2DenomVsCent->Write();
  v2Q1MidVsCent->Write();
  v2Q2MidVsCent->Write();
  v2EleNumVsCent->Write();
  v2EleDenomVsCent->Write();
  v2EleQ1MidVsCent->Write();
  v2EleQ2MidVsCent->Write();
    
  output->Close();

  timer.Stop();
  timer.Report();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: Z_EE_Channel <fileList> <job #> <total number of jobs>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);

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
   
  doZ2EE(listOfFiles, job);
  return 0; 
}

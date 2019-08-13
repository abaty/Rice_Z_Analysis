#include "include/centralityBin.h"
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

void doZ2EE(std::vector< std::string > files, int jobNumber, bool isTest){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  TH1::SetDefaultSumw2();
  ElectronEnergyScale energyScale = ElectronEnergyScale("MC");
  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronTnP eTnP = ElectronTnP();  

  MCReweight vzRW = MCReweight("resources/vzReweight.root");
  
  Settings s = Settings();

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();
  
  TRandom3 * r = new TRandom3();
  
  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  
  TH1D * recoEff_pt_pass[nBins];
  TH1D * recoEff_pt_net[nBins];
  TH1D * recoEff_pt[nBins];  
  TEfficiency * eff_pt[nBins];
  
  TH1D * recoEff_y_pass[nBins];
  TH1D * recoEff_y_net[nBins];
  TH1D * recoEff_y[nBins];  
  TEfficiency * eff_y[nBins];
  
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


  for(int k = 0; k<nBins; k++){
    massPeakOS[k] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[k] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_U_pass[k] = new TH2D(Form("recoEff_U_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_D_pass[k] = new TH2D(Form("recoEff_D_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_photonU_pass[k] = new TH2D(Form("recoEff_photonU_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_photonD_pass[k] = new TH2D(Form("recoEff_photonD_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_noSF_pass[k] = new TH2D(Form("recoEff_noSF_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_noSF_net[k] = new TH2D(Form("recoEff_noSF_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_pt_pass[k] = new TH1D(Form("recoEff_pt_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_net[k] = new TH1D(Form("recoEff_pt_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_y_pass[k] = new TH1D(Form("recoEff_y_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_y_net[k] = new TH1D(Form("recoEff_y_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_phi_pass[k] = new TH1D(Form("recoEff_phi_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_phi_net[k] = new TH1D(Form("recoEff_phi_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_cent_pass[k] = new TH1D(Form("recoEff_cent_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
    recoEff_cent_net[k] = new TH1D(Form("recoEff_cent_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
  
    
    recoEff_pt_pass_forReso_Reco[k] = new TH1D(Form("recoEff_pt_pass_forReso_Reco_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_RecoSmeared[k] = new TH1D(Form("recoEff_pt_pass_forReso_RecoSmeared_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_Gen[k] = new TH1D(Form("recoEff_pt_pass_forReso_Gen_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_ptReso[k] = new TH1D(Form("recoEff_ptReso_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.2,0.2);
    
    yReso[k] = new TH1D(Form("yReso_%d_%d",c.getCentBinLow(k), c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.1,0.1);
  }
  yReco = new TH1D("yReco",";y_{reco}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  yGen = new TH1D("yGen",";y_{gen}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
  yResponse = new TH2D("yResponse","y_{reco};y_{gen}",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);

  int nEle;
  int hiBin;
  float hiHF;
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
  if(isTest) nFiles = 5;
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
      double eventWeight = vzRW.reweightFactor( vz ) * c.findNcoll( hiBin );
      
      timer.StartSplit("Loading GEN electron tree");
      eTreeMC->GetEntry(i);

      timer.StartSplit("Gen Loop (PreElectronConfirmation)");
      int nGenElectronsFound = 0;
      
      TLorentzVector mom = TLorentzVector();
      bool foundGen = false;
      for(int j = 0; j<nMC; j++){
        //break out if you find a Z->tautau
        if( mcGMomPID->at(j) == 23 && TMath::Abs(mcMomPID->at(j)) == 15) break;
        
        //only look for electrons coming directly from a Z
        if( TMath::Abs(mcMomPID->at(j)) != 23) continue;
        
        //if it decays to muons, break out
        if( TMath::Abs(mcPID->at(j)) == 13) break;       
   
        //make sure it's an electron
        if( TMath::Abs(mcPID->at(j)) != 11) continue;       
        nGenElectronsFound++;

        //if they are not in our acceptance, break out
        if( mcPt->at(j) < s.minElectronPt ) break; 
        if( TMath::Abs(mcEta->at(j)) > s.maxZRapEle ) break; 
      
        timer.StartSplit("Gen Loop (PostElectronConfirmation)");

        //we found both daughters and they are in our acceptance, lets fill our histogram
        if( nGenElectronsFound == 2 ){
          //Fill denominator
          foundGen = true; 
          mom.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              recoEff_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_noSF_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_pt_net[k]->Fill( mom.Pt(), eventWeight);
              recoEff_y_net[k]->Fill( mom.Rapidity(), eventWeight);
              recoEff_cent_net[k]->Fill( hiBin/2.0, eventWeight);
              recoEff_phi_net[k]->Fill( mom.Phi(), eventWeight);

              //for the last few pt bins, halve the y binning so we have better stats
              if(mcMomPt->at(j) > s.zPtBins[s.nZPtBins - s.nPtBinsToRebinRapEff]){
                int bin = recoEff_net[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                if( bin%2 ==1){
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                } else {
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
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
          if( Zcand.Pt() < s.minPtCutForPhotons && acoplanarity < s.acoCutForPhotons ) passesAco[0] = false;
          if( Zcand.Pt() < s.minPtCutForPhotonsU && acoplanarity < s.acoCutForPhotonsU ) passesAco[1] = false;
          if( Zcand.Pt() < s.minPtCutForPhotonsD && acoplanarity < s.acoCutForPhotonsD ) passesAco[2] = false;
          
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
              if( mom.Pt() < s.zPtBins[ s.nZPtBins-1 ] && TMath::Abs( mom.Rapidity() ) < s.maxZRap ){
                if(passesAco[0]){

                  //90%
                  if(k==25){
                    yReco->Fill( Zcand.Rapidity(), eventWeight);
                    yGen->Fill( mom.Rapidity(), eventWeight);
                    yResponse->Fill(Zcand.Rapidity(), mom.Rapidity(), eventWeight);
                  }

                  massPeakOS[k]->Fill(Zcand.M(), eventWeight*scaleFactor);
                  recoEff_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                  recoEff_U_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactorU );
                  recoEff_D_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactorD );
                  recoEff_pt_pass[k]->Fill( mom.Pt(), eventWeight * scaleFactor);
                  recoEff_y_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                  recoEff_cent_pass[k]->Fill( hiBin/2.0, eventWeight * scaleFactor);
                  recoEff_phi_pass[k]->Fill( mom.Phi(), eventWeight * scaleFactor);
              
                  recoEff_pt_pass_forReso_Reco[k]->Fill(Zcand.Pt(), eventWeight);         
                  recoEff_pt_pass_forReso_RecoSmeared[k]->Fill(Zcand.Pt() * r->Gaus(1,0.05), eventWeight);         
                  recoEff_pt_pass_forReso_Gen[k]->Fill(mom.Pt(), eventWeight);         
                  recoEff_ptReso[k]->Fill( (Zcand.Pt() - mom.Pt()) / mom.Pt(), eventWeight);         
                
                  yReso[k]->Fill( Zcand.Rapidity() - mom.Rapidity() , eventWeight);            
          
                  //for the last few pt bins, halve the y binning so we have better stats
                  if(mom.Pt()> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
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
                  if(mom.Pt()> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
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
                  if(mom.Pt()> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
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
    forceConsistency(recoEff_y_pass[i], recoEff_y_net[i]);
    recoEff_y[i] = (TH1D*)recoEff_y_pass[i]->Clone(Form("recoEff_y_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_y[i]->Divide(recoEff_y_net[i]);
    recoEff_y[i]->SetDirectory(0);

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
    
    recoEff_y[i]->Write();
    recoEff_y_pass[i]->Write();
    recoEff_y_net[i]->Write();
    eff_y[i]->Write();
    
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
    massPeakSS[i]->Write();
    
    yReso[i]->Write();
  }
   
  yReco->Write();
  yGen->Write();
  yResponse->Write();
  yRatio = (TH1D*)yReco->Clone("yRatio");
  yRatio->Divide(yGen);
  yRatio->Write();
 
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

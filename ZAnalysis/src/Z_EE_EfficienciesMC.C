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

void doZ2EE(std::vector< std::string > files, int jobNumber){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  TH1::SetDefaultSumw2();
  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronTnP eTnP = ElectronTnP();  

  MCReweight vzRW = MCReweight("resources/vzReweight.root");
  
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoEff[nBins];
  TEfficiency * eff[nBins];

  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
  }

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
  std::vector< float > * eleHoverE = 0;
  std::vector< float > * eleD0 = 0;
  std::vector< float > * eleDz = 0;
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

  for(unsigned int f = 0; f<files.size(); f++){
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
    eTree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
    eTree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
    eTree->SetBranchAddress("eleSCEta",&eleSCEta);
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverE",&eleHoverE);
    eTree->SetBranchAddress("eleD0",&eleD0);
    eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("vz",&vz);
  
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
      
      double eventWeight = vzRW.reweightFactor( vz ) * c.findNcoll( hiBin );
      
      timer.StartSplit("Loading GEN electron tree");
      eTreeMC->GetEntry(i);

      timer.StartSplit("Gen Loop");
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

        //we found both daughters and they are in our acceptance, lets fill our histogram
        if( nGenElectronsFound == 2 ){
          //Fill denominator
          foundGen = true; 
          mom.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              recoEff_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );

              //for the last few pt bins, halve the y binning so we have better stats
              if(mcMomPt->at(j) > s.zPtBins[s.nZPtBins - s.nPtBinsToRebinRapEff]){
                int bin = recoEff_net[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                if( bin%2 ==1) recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                else           recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
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
        if(elePt->at(j)< s.minElectronPt) continue;
        if(TMath::Abs(eleSCEta->at(j)) > s.maxZRapEle) continue;
        //veto on transition region
        if(TMath::Abs(eleSCEta->at(j)) > 1.442 && TMath::Abs(eleSCEta->at(j)) < 1.556 ) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        float dEta = TMath::Abs( eledEtaAtVtx->at(j) );
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverE->at(j), eleEoverPInv->at(j), eleD0->at(j), eleDz->at(j))) continue;

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
          bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
          bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
          if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
  

          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " " << Zcand.Eta() << " " << Zcand.Phi() << " " << mom.Pt() << " " << mom.Eta() << " " << mom.Phi() <<  " isOS? " << (int)isOppositeSign << std::endl;
          if(!isOppositeSign) continue;
          if(moreThan2 && TMath::ACos(TMath::Cos(Zcand.Phi() - mom.Phi())) > 0.1 ) continue;
          //Looks like this is a good candidate match! Let's get the scale factor
          float scaleFactor = eTnP.getZSF(hiBin, elePt->at(goodElectrons.at(j)), eleSCEta->at(goodElectrons.at(j)), elePt->at(goodElectrons.at(j2)), eleSCEta->at(goodElectrons.at(j2)), 0) ;

          //Fill numerator (and apply the scale factor here!)
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
              if( mom.Pt() < s.zPtBins[ s.nZPtBins-1 ] && TMath::Abs( mom.Rapidity() ) < s.maxZRap ){
                recoEff_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
          
                //for the last few pt bins, halve the y binning so we have better stats
                if(mom.Pt()> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
                  int bin = recoEff_pass[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                  if( bin%2 ==1) recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactor );
                  else           recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactor );
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

  TFile * output = new TFile(Form("resources/Z2ee_EfficiencyMC_%d.root",jobNumber),"recreate");
  for(int i = 0; i<nBins; i++){
    recoEff_net[i]->Write();
    recoEff_pass[i]->Write();
    recoEff[i]->Write();
    eff[i]->Write();
  }
    
  output->Close();

  timer.Stop();
  timer.Report();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: Z_EE_EfficiencyMC <fileList> <job #> <total number of jobs>" << std::endl;
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

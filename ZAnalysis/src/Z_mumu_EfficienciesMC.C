#include "include/centralityBin.h"
#include "include/forceConsistency.h"
#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/MCReweight.h"
#include "include/MuonTnP.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TEfficiency.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

float getPhi(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2){
  TLorentzVector v1 = TLorentzVector();
  v1.SetPtEtaPhiM(pt1, eta1, phi1, 0.1056);
  TLorentzVector v2 = TLorentzVector();
  v2.SetPtEtaPhiM(pt2, eta2, phi2, 0.1056);

  return (v1+v2).Phi();
}

void doZ2mumuMC(std::vector< std::string > files){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  MuonTnP tnp = MuonTnP();
  MCReweight vzRW = MCReweight("resources/vzReweight.root");

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  int hiBin;

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
  
  TH2D * recoEff_noSF_pass[nBins];
  TH2D * recoEff_noSF_net[nBins];
  TH2D * recoEff_noSF[nBins];
  TEfficiency * eff_noSF[nBins];

  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    
    recoEff_noSF_pass[k] = new TH2D(Form("recoEff_noSF_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    recoEff_noSF_net[k] = new TH2D(Form("recoEff_noSF_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    
    recoEff_pt_pass[k] = new TH1D(Form("recoEff_pt_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_net[k] = new TH1D(Form("recoEff_pt_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_y_pass[k] = new TH1D(Form("recoEff_y_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_y_net[k] = new TH1D(Form("recoEff_y_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_phi_pass[k] = new TH1D(Form("recoEff_phi_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_phi_net[k] = new TH1D(Form("recoEff_phi_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_cent_pass[k] = new TH1D(Form("recoEff_cent_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
    recoEff_cent_net[k] = new TH1D(Form("recoEff_cent_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
  }

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),"dimucontana_mc"); 
    //for(unsigned int i = 0; i<v.GetEntries(); i++){
    for(unsigned int i = 0; i<100000; i++){
      v.GetEntry(i);
      hiBin = cb.getHiBinFromhiHFSides(v.HFsumETPlus() , v.HFsumETMinus() ,3);   
   
      if(i%1000==0) std::cout << i << "/" << v.GetEntries() << std::endl;
      
      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      double eventWeight = vzRW.reweightFactor( v.bestvtxZ() ) * c.findNcoll( hiBin );

      for(unsigned int j = 0; j<v.candSize_gen(); j++){

        //only look at gen Z's
        if(v.PID_gen()[j] != 23) continue;

        //require both legs to be in acceptance
        if( TMath::Abs( v.EtaD1_gen()[j] ) > s.maxZRap ) continue;
        if( TMath::Abs( v.EtaD2_gen()[j] ) > s.maxZRap ) continue;

        //require both muons to be > 20 GeV
        if( v.pTD1_gen()[j] < s.minMuonPt ) continue;
        if( v.pTD2_gen()[j] < s.minMuonPt ) continue;

        //Fill denominator 
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin( hiBin ,k)){
            recoEff_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );
            recoEff_noSF_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );
            recoEff_pt_net[k]->Fill( v.pT_gen()[j], eventWeight);
            recoEff_y_net[k]->Fill( v.y_gen()[j], eventWeight);
            recoEff_cent_net[k]->Fill( hiBin , eventWeight);
            recoEff_phi_net[k]->Fill( getPhi( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], v.PhiD2_gen()[j] ), eventWeight);

            //for the last few pt bins, halve the y binning so we have better stats
            if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins - s.nPtBinsToRebinRapEff]){
              int bin = recoEff_net[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
              if( bin%2 ==1){
                recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight );
                recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight );
              }else{
                recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight );
                recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight );
              }
            }
          }
        }

        //see if we have a matched candidate
        if(v.RecIdx_gen()[j] < 0 ) continue;
        int indx = v.RecIdx_gen()[j];

        //see if it passes all our selections
        if( v.mass()[indx] < s.zMassRange[0] || v.mass()[indx] > s.zMassRange[1]) continue; 
        if( !(v.pTD1()[indx] > s.minMuonPt )) continue;
        if( !(v.pTD2()[indx] > s.minMuonPt )) continue;
        if( TMath::Abs(v.EtaD1()[indx]) > s.maxZRap ) continue;
        if( TMath::Abs(v.EtaD2()[indx]) > s.maxZRap ) continue;
        if( !(v.tightCand(indx,"POG"))) continue;//tight Muon 1 && tight Muon 2      
        if( !(v.VtxProb()[indx] >0.001)) continue; 
 
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[indx] != v.chargeD2()[indx];
        if( !isOppositeSign ) continue;

        //Looks like this is a good candidate match! Let's get the scale factor
        float scaleFactor = tnp.getZSF( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], 0);

        //Fill numerator (and apply the scale factor here!)
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin( hiBin ,k)){
            //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
            if( v.pT_gen()[j] < s.zPtBins[ s.nZPtBins-1 ] && TMath::Abs( v.y_gen()[j] ) < s.maxZRap ){
              recoEff_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
              recoEff_noSF_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight  );
              recoEff_pt_pass[k]->Fill(v.pT_gen()[j], eventWeight * scaleFactor);         
              recoEff_y_pass[k]->Fill(v.y_gen()[j], eventWeight * scaleFactor);         
              recoEff_cent_pass[k]->Fill( hiBin /2.0, eventWeight * scaleFactor);         
              recoEff_phi_pass[k]->Fill(getPhi( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], v.PhiD2_gen()[j] ), eventWeight * scaleFactor);         
 
              //for the last few pt bins, halve the y binning so we have better stats
              if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
                int bin = recoEff_pass[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
                if( bin%2 ==1){
                  recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactor );
                  recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight );
                }else{
                  recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactor );
                  recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight  );
                }
              }
            }
          }
        }
      }
    }
  }

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

    recoEff_pass[i]->SetDirectory(0);
    recoEff_net[i]->SetDirectory(0);
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

    recoEff_noSF_pass[i]->SetDirectory(0);
    recoEff_noSF_net[i]->SetDirectory(0);
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

  TFile * output = new TFile("resources/Z2mumu_Efficiencies.root","recreate");
  for(int i = 0; i<nBins; i++){
    recoEff[i]->Write();
    recoEff_pass[i]->Write();
    recoEff_net[i]->Write();
    if(isConsistent.at(i)) eff[i]->Write();
    
    recoEff_noSF[i]->Write();
    recoEff_noSF_pass[i]->Write();
    recoEff_noSF_net[i]->Write();
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
  }

  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_mumu_EfficienciesMC <fileList>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

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
   
  doZ2mumuMC(listOfFiles);
  return 0; 
}

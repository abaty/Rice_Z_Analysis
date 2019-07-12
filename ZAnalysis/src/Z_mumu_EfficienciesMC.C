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


void doZ2mumuMC(std::vector< std::string > files){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  MuonTnP tnp = MuonTnP();
  MCReweight vzRW = MCReweight("resources/vzReweight.root");

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoEff[nBins];
  TEfficiency * eff[nBins];

  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
  }

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),"dimucontana_mc"); 
    for(unsigned int i = 0; i<v.GetEntries(); i++){
    //for(unsigned int i = 0; i<1000000; i++){
      v.GetEntry(i);
      
      if(i%1000==0) std::cout << i << "/" << v.GetEntries() << std::endl;
      
      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      double eventWeight = vzRW.reweightFactor( v.bestvtxZ() ) * c.findNcoll( v.centrality() );

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
          if(c.isInsideBin(v.centrality(),k)){
            recoEff_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );

            //for the last few pt bins, halve the y binning so we have better stats
            if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins - s.nPtBinsToRebinRapEff]){
              int bin = recoEff_net[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
              if( bin%2 ==1) recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight );
              else           recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight );
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
          if(c.isInsideBin(v.centrality(),k)){
            //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
            if( v.pT_gen()[j] < s.zPtBins[ s.nZPtBins-1 ] && TMath::Abs( v.y_gen()[j] ) < s.maxZRap ){
              recoEff_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
          
              //for the last few pt bins, halve the y binning so we have better stats
              if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
                int bin = recoEff_pass[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
                if( bin%2 ==1) recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactor );
                else           recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactor );
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

  TFile * output = new TFile("resources/Z2mumu_Efficiencies.root","recreate");
  for(int i = 0; i<nBins; i++){
    recoEff[i]->Write();
    recoEff_pass[i]->Write();
    recoEff_net[i]->Write();
    if(isConsistent.at(i)) eff[i]->Write();
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

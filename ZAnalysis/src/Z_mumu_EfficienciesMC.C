#include "include/centralityBin.h"
#include "include/forceConsistency.h"
#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/MCReweight.h"
#include "include/MuonTnP.h"
#include "include/MCWeightHelper.h"

//ROOT stuff
#include "TRandom3.h"
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

void doZ2mumuMC(std::vector< std::string > files, bool isTest){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  MuonTnP tnp = MuonTnP();
  MCReweight vzRW = MCReweight("resources/vzReweight.root","resources/centralityFlatteningWeight.root");

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  MCWeightHelper weightHelper = MCWeightHelper();

  TRandom3 * r = new TRandom3();

  int hiBin;

  TH1D * accept24_pt_pass[120];
  TH1D * accept24_y_pass[120];
  TH1D * accept24_yields_pass[120];
  TH1D * accept24_pt_net[120];
  TH1D * accept24_y_net[120];
  TH1D * accept24_yields_net[120];
  TH1D * accept21_pt_pass[120];
  TH1D * accept21_y_pass[120];
  TH1D * accept21_yields_pass[120];
  TH1D * accept21_pt_net[120];
  TH1D * accept21_y_net[120];
  TH1D * accept21_yields_net[120];

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
  
  TH1D * recoEff_pt_pass_forReso_Reco21[nBins];
  TH1D * recoEff_pt_pass_forReso_RecoSmeared21[nBins];
  TH1D * recoEff_pt_pass_forReso_Gen21[nBins];
  TH1D * recoEff_ptReso21[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_Reco21[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_RecoSmeared21[nBins];
  TH1D * recoEff_pt_pass_forReso_Ratio_NominalToSmeared21[nBins];

  TH1D * yReso[nBins];

  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_U_pass[k] = new TH2D(Form("recoEff_U_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_D_pass[k] = new TH2D(Form("recoEff_D_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_photonU_pass[k] = new TH2D(Form("recoEff_photonU_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_photonD_pass[k] = new TH2D(Form("recoEff_photonD_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_noSF_pass[k] = new TH2D(Form("recoEff_noSF_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_noSF_net[k] = new TH2D(Form("recoEff_noSF_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins4Eff-1,s.zPtBins4Eff);
    
    recoEff_pt_pass[k] = new TH1D(Form("recoEff_pt_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
    recoEff_pt_net[k] = new TH1D(Form("recoEff_pt_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins4Eff-1,s.zPtBins4Eff);
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
    
    recoEff_pt_pass_forReso_Reco21[k] = new TH1D(Form("recoEff_pt_pass_forReso_Reco21_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_RecoSmeared21[k] = new TH1D(Form("recoEff_pt_pass_forReso_RecoSmeared21_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_pass_forReso_Gen21[k] = new TH1D(Form("recoEff_pt_pass_forReso_Gen21_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_ptReso21[k] = new TH1D(Form("recoEff_ptReso21_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.2,0.2);

    yReso[k] = new TH1D(Form("yReso_%d_%d",c.getCentBinLow(k), c.getCentBinHigh(k)),";#frac{p_{T}^{reco}-p_{T}^{gen}}{p_{T}^{gen}}",40,-0.1,0.1);
  }
  for(int i = 0; i<weightHelper.getSize(); i++){
    accept24_pt_pass[i] = new TH1D(Form("accept24_pt_pass_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept24_pt_net[i] = new TH1D(Form("accept24_pt_net_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept24_y_pass[i] = new TH1D(Form("accept24_y_pass_%d",i),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    accept24_y_net[i] = new TH1D(Form("accept24_y_net_%d",i),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    accept24_yields_pass[i] = new TH1D(Form("accept24_yields_pass_%d",i),"",1,0,1);
    accept24_yields_net[i] = new TH1D(Form("accept24_yields_net_%d",i),"",1,0,1);
    accept21_pt_pass[i] = new TH1D(Form("accept21_pt_pass_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept21_pt_net[i] = new TH1D(Form("accept21_pt_net_%d",i),"",s.nZPtBins-1,s.zPtBins);
    accept21_y_pass[i] = new TH1D(Form("accept21_y_pass_%d",i),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
    accept21_y_net[i] = new TH1D(Form("accept21_y_net_%d",i),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle);
    accept21_yields_pass[i] = new TH1D(Form("accept21_yields_pass_%d",i),"",1,0,1);
    accept21_yields_net[i] = new TH1D(Form("accept21_yields_net_%d",i),"",1,0,1);
  }

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),"dimucontana_mc"); 
    unsigned int maxEvt = isTest ? 100000 : v.GetEntries();
    for(unsigned int i = 0; i<maxEvt; i++){
      v.GetEntry(i);
   
      if(i%1000==0) std::cout << i << "/" << v.GetEntries() << std::endl;
      
      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;
      hiBin = cb.getHiBinFromhiHFSides(v.HFsumETPlus() , v.HFsumETMinus() ,3);   

      double acceptWeight[120];
      for(int w = 0; w<weightHelper.getSize(); w++){
        acceptWeight[w] = vzRW.reweightFactor( v.bestvtxZ() ) * (v.weightLHE_gen(weightHelper.getIndx(w))/10000.0);
      }
      double eventWeight = vzRW.reweightFactor( v.bestvtxZ() ) * vzRW.reweightFactorCent(hiBin) * c.findNcoll( hiBin ) * (v.weightLHE_gen(1080)/10000.0);//1080 is EPPS16

      for(unsigned int j = 0; j<v.candSize_gen(); j++){

        //only look at gen Z's
        if(v.PID_gen()[j] != 23) continue;
        if(v.DecayID_gen()[j] != 23) continue;

        //rapidity cut and acceptance stuff
        if(TMath::Abs( v.y_gen()[j] ) > s.maxZRap) continue;

        for(int w = 0; w<weightHelper.getSize(); w++){
          accept24_yields_net[w]->Fill( 0.5, acceptWeight[w]); 
          accept24_y_net[w]->Fill( v.y_gen()[j], acceptWeight[w]); 
          accept24_pt_net[w]->Fill( v.pT_gen()[j], acceptWeight[w]); 
        }
        
        if(TMath::Abs( v.y_gen()[j] ) < s.maxZRapEle){
          for(int w = 0; w<weightHelper.getSize(); w++){
            accept21_yields_net[w]->Fill( 0.5, acceptWeight[w]); 
            accept21_y_net[w]->Fill( v.y_gen()[j], acceptWeight[w]); 
            accept21_pt_net[w]->Fill( v.pT_gen()[j], acceptWeight[w]); 
          }
        }
 
        //require both legs to be in acceptance
        if( TMath::Abs( v.EtaD1_gen()[j] ) > s.maxZRap ) continue;
        if( TMath::Abs( v.EtaD2_gen()[j] ) > s.maxZRap ) continue;

        //require both muons to be > 20 GeV
        if( v.pTD1_gen()[j] < s.minMuonPt ) continue;
        if( v.pTD2_gen()[j] < s.minMuonPt ) continue;

        for(int w = 0; w<weightHelper.getSize(); w++){
          accept24_yields_pass[w]->Fill(0.5, acceptWeight[w]);
          accept24_y_pass[w]->Fill(v.y_gen()[j], acceptWeight[w]);
          accept24_pt_pass[w]->Fill(v.pT_gen()[j], acceptWeight[w]);
        }

        if( TMath::Abs( v.y_gen()[j] ) < s.maxZRapEle && TMath::Abs( v.EtaD1_gen()[j] ) < s.maxZRapEle && TMath::Abs( v.EtaD2_gen()[j] ) < s.maxZRapEle){
          for(int w = 0; w<weightHelper.getSize(); w++){
            accept21_yields_pass[w]->Fill(0.5, acceptWeight[w]);
            accept21_y_pass[w]->Fill(v.y_gen()[j], acceptWeight[w]);
            accept21_pt_pass[w]->Fill(v.pT_gen()[j], acceptWeight[w]);
          }
        }

        //Fill denominator 
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin( hiBin ,k)){
            recoEff_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );
            recoEff_noSF_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );
            recoEff_pt_net[k]->Fill( v.pT_gen()[j], eventWeight);
            recoEff_y_net[k]->Fill( v.y_gen()[j], eventWeight);
            recoEff_cent_net[k]->Fill( hiBin/2.0 , eventWeight);
            recoEff_phi_net[k]->Fill( getPhi( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], v.PhiD2_gen()[j] ), eventWeight);

            //for the last few pt bins, halve the y binning so we have better stats
            if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins4Eff - s.nPtBinsToRebinRapEff]){
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

        //we have 3 different aco cuts for variations
        bool passesAco[3] = {1 , 1, 1};
        float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( v.PhiD1()[indx] - v.PhiD2()[indx] )))/TMath::Pi(); 
        if(v.pT()[indx] < s.minPtCutForPhotons && acoplanarity < s.acoCutForPhotons) passesAco[0] = false;
        if(v.pT()[indx] < s.minPtCutForPhotonsU && acoplanarity < s.acoCutForPhotonsU) passesAco[1] = false;
        if(v.pT()[indx] < s.minPtCutForPhotonsD && acoplanarity < s.acoCutForPhotonsD) passesAco[2] = false;
 
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[indx] != v.chargeD2()[indx];
        if( !isOppositeSign ) continue;

        //Looks like this is a good candidate match! Let's get the scale factor
        float scaleFactor = tnp.getZSF( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], 0);
        float scaleFactorU = tnp.getZSF( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], 1);
        float scaleFactorD = tnp.getZSF( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], -1);

        //Fill numerator (and apply the scale factor here!)
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin( hiBin ,k)){
            //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
            if( v.pT_gen()[j] < s.zPtBins[ s.nZPtBins4Eff-1 ] && TMath::Abs( v.y_gen()[j] ) < s.maxZRap ){
              if(passesAco[0]){
                recoEff_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
                recoEff_U_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactorU );
                recoEff_D_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactorD );
                recoEff_noSF_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight  );
                recoEff_pt_pass[k]->Fill(v.pT_gen()[j], eventWeight * scaleFactor);         
                recoEff_y_pass[k]->Fill(v.y_gen()[j], eventWeight * scaleFactor);         
                recoEff_cent_pass[k]->Fill( hiBin /2.0, eventWeight * scaleFactor);         
                recoEff_phi_pass[k]->Fill(getPhi( v.pTD1_gen()[j], v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD2_gen()[j], v.EtaD2_gen()[j], v.PhiD2_gen()[j] ), eventWeight * scaleFactor);         
                
                recoEff_pt_pass_forReso_Reco[k]->Fill(v.pT()[indx], eventWeight);         
                recoEff_pt_pass_forReso_RecoSmeared[k]->Fill(v.pT()[indx] * r->Gaus(1,0.05), eventWeight);         
                recoEff_pt_pass_forReso_Gen[k]->Fill(v.pT_gen()[j], eventWeight);         
                recoEff_ptReso[k]->Fill( (v.pT()[indx] - v.pT_gen()[j]) / v.pT_gen()[j], eventWeight);         
                yReso[k]->Fill( v.y()[indx] - v.y_gen()[j] , eventWeight);            
 
                if(TMath::Abs(v.EtaD1()[indx]) < s.maxZRapEle && TMath::Abs(v.EtaD2()[indx]) < s.maxZRapEle){   
                  recoEff_pt_pass_forReso_Reco21[k]->Fill(v.pT()[indx], eventWeight);         
                  recoEff_pt_pass_forReso_RecoSmeared21[k]->Fill(v.pT()[indx] * r->Gaus(1,0.05), eventWeight);         
                  recoEff_pt_pass_forReso_Gen21[k]->Fill(v.pT_gen()[j], eventWeight);         
                  recoEff_ptReso21[k]->Fill( (v.pT()[indx] - v.pT_gen()[j]) / v.pT_gen()[j], eventWeight); 
                }        
   
                //for the last few pt bins, halve the y binning so we have better stats
                if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                  int bin = recoEff_pass[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
                  if( bin%2 ==1){
                    recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactor );
                    recoEff_U_pass[k]->Fill( recoEff_U_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactorU );
                    recoEff_D_pass[k]->Fill( recoEff_D_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactorD );
                    recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight );
                  }else{
                    recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactor );
                    recoEff_U_pass[k]->Fill( recoEff_U_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactorU );
                    recoEff_D_pass[k]->Fill( recoEff_D_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactorD );
                    recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight  );
                  }
                }
              }//passes aco if statement here
              if(passesAco[1]){
                recoEff_photonU_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
                if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                  int bin = recoEff_pass[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
                  if( bin%2 ==1){
                    recoEff_photonU_pass[k]->Fill( recoEff_photonU_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactor );
                  }else{
                    recoEff_photonU_pass[k]->Fill( recoEff_photonU_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactor );
                  }
                }
              }
              if(passesAco[2]){
                recoEff_photonD_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
                if(v.pT_gen()[j]> s.zPtBins[s.nZPtBins4Eff- s.nPtBinsToRebinRapEff]){
                  int bin = recoEff_pass[k]->GetXaxis()->FindBin(v.y_gen()[j]); 
                  if( bin%2 ==1){
                    recoEff_photonD_pass[k]->Fill( recoEff_photonD_pass[k]->GetXaxis()->GetBinCenter(bin+1), v.pT_gen()[j], eventWeight * scaleFactor );
                  }else{
                    recoEff_photonD_pass[k]->Fill( recoEff_photonD_pass[k]->GetXaxis()->GetBinCenter(bin-1), v.pT_gen()[j], eventWeight * scaleFactor );
                  }
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

  TFile * output;
  if(isTest) output = new TFile("resources/Z2mumu_Efficiencies_TEST.root","recreate");
  else       output = new TFile("resources/Z2mumu_Efficiencies.root","recreate");
  for(int i = 0; i<nBins; i++){
    recoEff[i]->Write();
    recoEff_pass[i]->Write();
    recoEff_net[i]->Write();
    if(isConsistent.at(i)) eff[i]->Write();
    
    recoEff_U[i]->Write();
    recoEff_U_pass[i]->Write();
    eff_U[i]->Write();
    
    recoEff_D[i]->Write();
    recoEff_D_pass[i]->Write();
    eff_D[i]->Write();
    
    recoEff_photonU[i]->Write();
    recoEff_photonU_pass[i]->Write();
    eff_photonU[i]->Write();
    
    recoEff_photonD[i]->Write();
    recoEff_photonD_pass[i]->Write();
    eff_photonD[i]->Write();
    
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
    
    recoEff_pt_pass_forReso_Ratio_Reco21[i] = (TH1D*) recoEff_pt_pass_forReso_Reco21[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_Reco21_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_Reco21[i]->Divide(recoEff_pt_pass_forReso_Gen21[i]);
    recoEff_pt_pass_forReso_Ratio_Reco21[i]->Write();
    recoEff_pt_pass_forReso_Ratio_RecoSmeared21[i] = (TH1D*) recoEff_pt_pass_forReso_RecoSmeared21[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_RecoSmeared21_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_RecoSmeared21[i]->Divide(recoEff_pt_pass_forReso_Gen21[i]);
    recoEff_pt_pass_forReso_Ratio_RecoSmeared21[i]->Write();

    recoEff_pt_pass_forReso_Ratio_NominalToSmeared21[i] = (TH1D*) recoEff_pt_pass_forReso_Ratio_RecoSmeared21[i]->Clone(Form("recoEff_pt_pass_forReso_Ratio_SmearedToNominal21_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt_pass_forReso_Ratio_NominalToSmeared21[i]->Divide(recoEff_pt_pass_forReso_Ratio_Reco21[i]);
    recoEff_pt_pass_forReso_Ratio_NominalToSmeared21[i]->Write();

    recoEff_pt_pass_forReso_Reco21[i]->Write();
    recoEff_pt_pass_forReso_RecoSmeared21[i]->Write();
    recoEff_pt_pass_forReso_Gen21[i]->Write();
    recoEff_ptReso21[i]->Write();

    yReso[i]->Write();
  }

  TH1D * accept24_yields_ratio[120];
  TH1D * accept24_y_ratio[120];
  TH1D * accept24_pt_ratio[120];
  TH1D * accept24_scaleVariation_yields_ratio[120];
  TH1D * accept24_scaleVariation_y_ratio[120];
  TH1D * accept24_scaleVariation_pt_ratio[120];
  TH1D * accept24_nPDFVariation_yields_ratio[120];
  TH1D * accept24_nPDFVariation_y_ratio[120];
  TH1D * accept24_nPDFVariation_pt_ratio[120];
  TH1D * accept24_nPDFVariationMax_yields_ratio;
  TH1D * accept24_nPDFVariationMax_y_ratio;
  TH1D * accept24_nPDFVariationMax_pt_ratio;
  for(int w = 0; w<weightHelper.getSize(); w++){
    accept24_yields_ratio[w] = (TH1D*)accept24_yields_pass[w]->Clone(Form("accept24_yields_ratio_%d",w));
    accept24_yields_ratio[w]->Divide(accept24_yields_net[w]);
    accept24_yields_ratio[w]->Write();
    
    accept24_y_ratio[w] = (TH1D*)accept24_y_pass[w]->Clone(Form("accept24_y_ratio_%d",w));
    accept24_y_ratio[w]->Divide(accept24_y_net[w]);
    accept24_y_ratio[w]->Write();
    
    accept24_pt_ratio[w] = (TH1D*)accept24_pt_pass[w]->Clone(Form("accept24_pt_ratio_%d",w));
    accept24_pt_ratio[w]->Divide(accept24_pt_net[w]);
    accept24_pt_ratio[w]->Write();

    //accept24_yields_pass[w]->Write();
    //accept24_yields_net[w]->Write();
    //accept24_pt_pass[w]->Write();
    //accept24_pt_net[w]->Write();
    //accept24_y_pass[w]->Write();
    //accept24_y_net[w]->Write();

    //scale variations
    if(w>=3 && w<=8){
      accept24_scaleVariation_yields_ratio[w] = (TH1D*) accept24_yields_ratio[w]->Clone(Form("accept24_scaleVariation_yields_%d",w));
      accept24_scaleVariation_yields_ratio[w]->Divide( accept24_yields_ratio[2] );
      accept24_scaleVariation_yields_ratio[w]->Write();
      accept24_scaleVariation_y_ratio[w] = (TH1D*) accept24_y_ratio[w]->Clone(Form("accept24_scaleVariation_y_%d",w));
      accept24_scaleVariation_y_ratio[w]->Divide( accept24_y_ratio[2] );
      accept24_scaleVariation_y_ratio[w]->Write();
      accept24_scaleVariation_pt_ratio[w] = (TH1D*) accept24_pt_ratio[w]->Clone(Form("accept24_scaleVariation_pt_%d",w));
      accept24_scaleVariation_pt_ratio[w]->Divide( accept24_pt_ratio[2] );
      accept24_scaleVariation_pt_ratio[w]->Write();
    }
    if(w==1 || w==2 || w>8){
      accept24_nPDFVariation_yields_ratio[w] = (TH1D*) accept24_yields_ratio[w]->Clone(Form("accept24_nPDFVariation_yields_%d",w));
      accept24_nPDFVariation_yields_ratio[w]->Divide( accept24_yields_ratio[0] );
      accept24_nPDFVariation_yields_ratio[w]->Write();
      accept24_nPDFVariation_y_ratio[w] = (TH1D*) accept24_y_ratio[w]->Clone(Form("accept24_nPDFVariation_y_%d",w));
      accept24_nPDFVariation_y_ratio[w]->Divide( accept24_y_ratio[0] );
      accept24_nPDFVariation_y_ratio[w]->Write();
      accept24_nPDFVariation_pt_ratio[w] = (TH1D*) accept24_pt_ratio[w]->Clone(Form("accept24_nPDFVariation_pt_%d",w));
      accept24_nPDFVariation_pt_ratio[w]->Divide( accept24_pt_ratio[0] );
      accept24_nPDFVariation_pt_ratio[w]->Write();
 
      if(w>8){
        if(w==9){
          accept24_nPDFVariationMax_yields_ratio = (TH1D*) accept24_nPDFVariation_yields_ratio[w]->Clone("accept24_nPDFVariationMax_yields_ratio");
          accept24_nPDFVariationMax_yields_ratio->Reset();
          accept24_nPDFVariationMax_y_ratio = (TH1D*) accept24_nPDFVariation_y_ratio[w]->Clone("accept24_nPDFVariationMax_y_ratio");
          accept24_nPDFVariationMax_y_ratio->Reset();
          accept24_nPDFVariationMax_pt_ratio = (TH1D*) accept24_nPDFVariation_pt_ratio[w]->Clone("accept24_nPDFVariationMax_pt_ratio");
          accept24_nPDFVariationMax_pt_ratio->Reset();
        }
        for(int bin = 0; bin<accept24_nPDFVariation_yields_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept24_nPDFVariation_yields_ratio[w]->GetBinContent(bin) );
          if( content > accept24_nPDFVariationMax_yields_ratio->GetBinContent(bin)){
            accept24_nPDFVariationMax_yields_ratio->SetBinContent(bin, content);
          }
        }
        for(int bin = 0; bin<accept24_nPDFVariation_y_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept24_nPDFVariation_y_ratio[w]->GetBinContent(bin) );
          if( content > accept24_nPDFVariationMax_y_ratio->GetBinContent(bin)){
            accept24_nPDFVariationMax_y_ratio->SetBinContent(bin, content);
          }
        }
        for(int bin = 0; bin<accept24_nPDFVariation_pt_ratio[w]->GetSize(); bin++){
          float content = TMath::Abs( 1 - accept24_nPDFVariation_pt_ratio[w]->GetBinContent(bin) );
          if( content > accept24_nPDFVariationMax_pt_ratio->GetBinContent(bin)){
            accept24_nPDFVariationMax_pt_ratio->SetBinContent(bin, content);
          }
        }
      }
    }
  }
  accept24_nPDFVariationMax_yields_ratio->Scale(1.0/1.645);//envelope gives 90%, scale down to 68% coverage
  accept24_nPDFVariationMax_y_ratio->Scale(1.0/1.645);
  accept24_nPDFVariationMax_pt_ratio->Scale(1.0/1.645);
  accept24_nPDFVariationMax_yields_ratio->Write();
  accept24_nPDFVariationMax_y_ratio->Write();
  accept24_nPDFVariationMax_pt_ratio->Write();

  TH1D * accept24_pdfTypeVariation_yields_ratio = (TH1D*) accept24_nPDFVariation_yields_ratio[1]->Clone("accept24_pdfTypeVariation_yields_ratio");
  TH1D * accept24_totalUncert_yields = (TH1D*) accept24_nPDFVariationMax_yields_ratio->Clone("accept24_totalUncert_yields");
  for(int bin = 0; bin<accept24_pdfTypeVariation_yields_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept24_pdfTypeVariation_yields_ratio->GetBinContent(bin) );
    accept24_pdfTypeVariation_yields_ratio->SetBinContent(bin, content);

    float content2 = accept24_totalUncert_yields->GetBinContent(bin);
    accept24_totalUncert_yields->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept24_pdfTypeVariation_yields_ratio->Write();
  accept24_totalUncert_yields->Write();
  
  TH1D * accept24_pdfTypeVariation_y_ratio = (TH1D*) accept24_nPDFVariation_y_ratio[1]->Clone("accept24_pdfTypeVariation_y_ratio");
  TH1D * accept24_totalUncert_y = (TH1D*) accept24_nPDFVariationMax_y_ratio->Clone("accept24_totalUncert_y");
  for(int bin = 0; bin<accept24_pdfTypeVariation_y_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept24_pdfTypeVariation_y_ratio->GetBinContent(bin) );
    accept24_pdfTypeVariation_y_ratio->SetBinContent(bin, content);
    
    float content2 = accept24_totalUncert_y->GetBinContent(bin);
    accept24_totalUncert_y->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept24_pdfTypeVariation_y_ratio->Write();
  accept24_totalUncert_y->Write();  

  TH1D * accept24_pdfTypeVariation_pt_ratio = (TH1D*) accept24_nPDFVariation_pt_ratio[1]->Clone("accept24_pdfTypeVariation_pt_ratio");
  TH1D * accept24_totalUncert_pt = (TH1D*) accept24_nPDFVariationMax_pt_ratio->Clone("accept24_totalUncert_pt");
  for(int bin = 0; bin<accept24_pdfTypeVariation_pt_ratio->GetSize(); bin++){
    float content = TMath::Abs( 1 - accept24_pdfTypeVariation_pt_ratio->GetBinContent(bin) );
    accept24_pdfTypeVariation_pt_ratio->SetBinContent(bin, content);
    
    float content2 = accept24_totalUncert_pt->GetBinContent(bin);
    accept24_totalUncert_pt->SetBinContent(bin, TMath::Sqrt(content* content + content2 * content2));
  }
  accept24_pdfTypeVariation_pt_ratio->Write();
  accept24_totalUncert_pt->Write();

  
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

  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: Z_mumu_EfficienciesMC <fileList> <isTest>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  bool isTest = (bool)std::atoi(argv[2]);
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
   
  doZ2mumuMC(listOfFiles, isTest);
  return 0; 
}

#!/bin/bash
source environment.sh

./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBin1_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 21hiBin1
./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBin2_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 21hiBin2
./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBin0_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 21
#for ZDC crosscehck
#./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBinZDC0_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root 1 21hiBin1
#./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBinZDC0_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root 1 21hiBin2
#./bin/massPeakPlots_BkgSub.exe Z2mumu_21_Data_hiBinZDC0_job0.root Z2mumu_MC_21_DY_job0.root Z2mumu_MC_21_TTbar_job0.root Z2mumu_MC_21_WJet_job0.root 1 21

./bin/massPeakPlots_BkgSub.exe Z2mumu_24_Data_hiBin1_job0.root Z2mumu_MC_24_DY_job0.root Z2mumu_MC_24_TTbar_job0.root Z2mumu_MC_24_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 24hiBin1
./bin/massPeakPlots_BkgSub.exe Z2mumu_24_Data_hiBin2_job0.root Z2mumu_MC_24_DY_job0.root Z2mumu_MC_24_TTbar_job0.root Z2mumu_MC_24_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 24hiBin2
./bin/massPeakPlots_BkgSub.exe Z2mumu_24_Data_hiBin0_job0.root Z2mumu_MC_24_DY_job0.root Z2mumu_MC_24_TTbar_job0.root Z2mumu_MC_24_WJet_job0.root resources/Z2mumu_Efficiencies.root 1 24

./bin/massPeakPlots_BkgSub.exe Z2ee_hiBin1.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root resources/Z2ee_EfficiencyMC_0.root 0 eehiBin1
./bin/massPeakPlots_BkgSub.exe Z2ee_hiBin2.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root resources/Z2ee_EfficiencyMC_0.root 0 eehiBin2
./bin/massPeakPlots_BkgSub.exe Z2ee.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root resources/Z2ee_EfficiencyMC_0.root 0 ee
#./bin/massPeakPlots_BkgSub.exe Z2ee_Data_hiBinZDC0.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root 0 eehiBin1
#./bin/massPeakPlots_BkgSub.exe Z2ee_Data_hiBinZDC0.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root 0 eehiBin2
#./bin/massPeakPlots_BkgSub.exe Z2ee_Data_hiBinZDC0.root Z2ee_MC_DY.root Z2ee_MC_TTbar.root Z2ee_MC_WJet.root 0 ee


./bin/systematics.exe backgroundSubtraction_ee_isMu0.root backgroundSubtraction_eehiBin1_isMu0.root backgroundSubtraction_eehiBin2_isMu0.root resources/Z2ee_EfficiencyMC_0.root ee 0 1
./bin/systematics.exe backgroundSubtraction_21_isMu1.root backgroundSubtraction_21hiBin1_isMu1.root backgroundSubtraction_21hiBin2_isMu1.root resources/Z2mumu_Efficiencies.root mu21 1 0
./bin/systematics.exe backgroundSubtraction_24_isMu1.root backgroundSubtraction_24hiBin1_isMu1.root backgroundSubtraction_24hiBin2_isMu1.root resources/Z2mumu_Efficiencies.root mu24 0 0

./bin/systematicsV2.exe Z2ee.root Z2ee_hiBin1.root Z2ee_hiBin2.root 1
./bin/systematicsV2.exe Z2mumu_24_Data_hiBin0_job0.root Z2mumu_24_Data_hiBin1_job0.root Z2mumu_24_Data_hiBin2_job0.root 0

./bin/prettyPlots_yields.exe backgroundSubtraction_ee_isMu0.root backgroundSubtraction_21_isMu1.root backgroundSubtraction_24_isMu1.root systematics_ee_isMu210.root systematics_mu21_isMu211.root systematics_mu24_isMu210.root 
./bin/prettyPlots_pTy.exe backgroundSubtraction_ee_isMu0.root backgroundSubtraction_21_isMu1.root backgroundSubtraction_24_isMu1.root systematics_ee_isMu210.root systematics_mu21_isMu211.root systematics_mu24_isMu210.root 
./bin/v2Plots.exe Z2ee.root Z2mumu_24_Data_hiBin0_job0.root systematics_v2_isEE1.root  systematics_v2_isEE0.root 


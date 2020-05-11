#ifndef SETTINGS
#define SETTINGS


class Settings{
public:

  double sigmaMBPbPb = 7.644;

  double Nmb = 11536000857;
  double NmbUp = 1.012*11536000857;
  double NmbDown = 0.988*11536000857;

  double muLumi = 1682.8;
  double eLumi = 1670.7 ;
  double netLumi = 1695.6;


  double minMuonPt = 20;
  double minElectronPt = 20;

  int nZMassBins = 30;
  double zMassRange[2] = {60,120};

  static const int nZPtBins = 14;
  double zPtBins[nZPtBins] = {0,1.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,90.0,120.0,150.0,200.0};
  
  static const int nZPtBins4Eff = 9;
  double zPtBins4Eff[nZPtBins4Eff] = {0,1.0,3.0,5.0,10.0,20.0,40.0,70.0,200.0};
  int nPtBinsToRebinRapEff = 3;

  static const int nZRapBins = 16;
  double maxZRap = 2.4;
  static const int nZRapBinsEle = 14;
  double maxZRapEle = 2.1;

  int nMuEffBinsEta = 48;
  static const int nMuEffBinsPt = 11;
  double muEffBinsPt[nMuEffBinsPt+1] = {20,25,30,35,40,45,50,60,80,120,160,200};

  double minPtCutForPhotons = 1.25;
  double acoCutForPhotons = 0.001;
  double minPtCutForPhotonsU = 1.75;
  double acoCutForPhotonsU = 0.0015;
  double minPtCutForPhotonsD = 1.0;
  double acoCutForPhotonsD = 0.00075;
  
  double minPtCutForPhotons_ELE = 2.5;
  double acoCutForPhotons_ELE = 0.001;
  double minPtCutForPhotonsU_ELE = 3.0;
  double acoCutForPhotonsU_ELE = 0.0015;
  double minPtCutForPhotonsD_ELE = 2.0;
  double acoCutForPhotonsD_ELE = 0.00075;

  //MC xsections (/pb)
  double ttbar_XS = 69.0;
  double Wjet_XS = 21159;
  double DY_XS = 2010;

  
private:

};

#endif

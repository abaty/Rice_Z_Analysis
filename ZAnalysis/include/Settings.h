#ifndef SETTINGS
#define SETTINGS


class Settings{
public:

  double Nmb = 9.57519e+09;
  double NmbUp = 9.62460e+09;
  double NmbDown = 9.47799e+09;

  double muLumi = 1606.06;
  double eLumi = 1594.214;
  double netLumi = 1618.466;


  double minMuonPt = 20;
  double minElectronPt = 20;

  int nZMassBins = 30;
  double zMassRange[2] = {60,120};

  static const int nZPtBins = 14;
  double zPtBins[nZPtBins] = {0,1.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,90.0,120.0,150.0,200.0};
  int nPtBinsToRebinRapEff = 6;

  static const int nZRapBins = 16;
  double maxZRap = 2.4;
  static const int nZRapBinsEle = 14;
  double maxZRapEle = 2.1;

  int nMuEffBinsEta = 48;
  static const int nMuEffBinsPt = 11;
  double muEffBinsPt[nMuEffBinsPt+1] = {20,25,30,35,40,45,50,60,80,120,160,200};

  double minPtCutForPhotons = 0.5;
  double acoCutForPhotons = 0.001;
  
  double minPtCutForPhotonsU = 0.75;
  double acoCutForPhotonsU = 0.0015;
  
  double minPtCutForPhotonsD = 0.25;
  double acoCutForPhotonsD = 0.0005;

  //MC xsections (/pb)
  double ttbar_XS = 69.0;
  double Wjet_XS = 21159;
  double DY_XS = 2010;

  
private:

};

#endif

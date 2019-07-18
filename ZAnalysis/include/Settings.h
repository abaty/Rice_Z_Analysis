#ifndef SETTINGS
#define SETTINGS


class Settings{
public:

  double minMuonPt = 20;
  double minElectronPt = 20;

  int nZMassBins = 30;
  double zMassRange[2] = {60,120};

  static const int nZPtBins = 15;
  double zPtBins[nZPtBins] = {0,0.5,1.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,90.0,120.0,150.0,200.0};
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

  //MC xsections (/pb)
  double ttbar_XS = 69.0;
  double Wjet_XS = 21159;
  double DY_XS = 2010;

  
private:

};

#endif

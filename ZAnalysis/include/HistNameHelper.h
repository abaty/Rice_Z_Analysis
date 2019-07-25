#ifndef HISTHELPER
#define HISTHELPER

#include "TColor.h"
#include "TBox.h"
#include "TMath.h"
#include "TH1D.h"
#include <string>
#include <vector>

class HistNameHelper{
  public:

  HistNameHelper();
  void makeDifferential(TH1D * h);
  void subtractOneAndAbs(TH1D * h);
  void subtractOneAndAbsAndSymmetrize(TH1D * h, TH1D * h2);
  void absAndSymmetrize(TH1D * h, TH1D * h2);
  void addInQuadrature2(TH1D * h1, TH1D * h2);
  void addInQuadrature3(TH1D * h1, TH1D * h2, TH1D * h3);
  void addInQuadrature4(TH1D * h1, TH1D * h2, TH1D * h3, TH1D * h4);

  void drawBoxAbsolute(TH1D * h, int i, TBox * b, float absError, float width, Color_t color, bool doDraw = true);
  ~HistNameHelper();

  std::vector< std::string > name;
  std::vector< std::string > variationName;
};

HistNameHelper::HistNameHelper(){
  name.push_back("massPeak");
  name.push_back("pT");
  name.push_back("y");
  name.push_back("yield");


  variationName = {"","_var1","_var2","_var3","_var4"};
}

HistNameHelper::~HistNameHelper(){

}

void HistNameHelper::drawBoxAbsolute(TH1D * h, int i, TBox * b, float absError, float width, Color_t color, bool doDraw){
  b = new TBox(h->GetBinCenter(i)-width,h->GetBinContent(i)-absError,h->GetBinCenter(i)+width,h->GetBinContent(i)+absError);
  b->SetLineColor(color);
  b->SetFillStyle(0);
  if(doDraw) b->Draw("same");
}

void HistNameHelper::makeDifferential(TH1D* h){
  for(int j = 1; j<h->GetNbinsX()+1; j++){
    h->SetBinContent(j,h->GetBinContent(j)/h->GetBinWidth(j));
    h->SetBinError(j,h->GetBinError(j)/h->GetBinWidth(j));
  }
  return;
}

//takes a histogram, subtracts off 1 and then takes the abs of it
void HistNameHelper::subtractOneAndAbs(TH1D* h){
  for(int j = 1; j<h->GetNbinsX()+1; j++){
    h->SetBinContent(j,TMath::Abs(h->GetBinContent(j)-1));
    h->SetBinError(j,0);
  }
  return;
}

//sets the first histogram to the max of |h-1| for both histograms 
void HistNameHelper::subtractOneAndAbsAndSymmetrize(TH1D* h, TH1D * h2){
  for(int j = 1; j<h->GetNbinsX()+1; j++){
    h->SetBinContent(j, TMath::Max( TMath::Abs(h->GetBinContent(j)-1), TMath::Abs(h2->GetBinContent(j)-1) ));
    h->SetBinError(j,0);
  }
  return;
}

//sets the first histogram to the max of |h-1| for both histograms 
void HistNameHelper::absAndSymmetrize(TH1D* h, TH1D * h2){
  for(int j = 1; j<h->GetNbinsX()+1; j++){
    h->SetBinContent(j, TMath::Max( TMath::Abs(h->GetBinContent(j)), TMath::Abs(h2->GetBinContent(j)) ));
    h->SetBinError(j,0);
  }
  return;
}

void HistNameHelper::addInQuadrature2(TH1D *h1, TH1D *h2){
  for(int j = 1; j<h1->GetNbinsX()+1; j++){
    h1->SetBinContent(j, TMath::Sqrt(TMath::Power(h1->GetBinContent(j),2) + TMath::Power(h2->GetBinContent(j),2)) );
    h1->SetBinError(j,0);
  }
  return;

}
void HistNameHelper::addInQuadrature3(TH1D *h1, TH1D *h2, TH1D *h3){
  for(int j = 1; j<h1->GetNbinsX()+1; j++){
    h1->SetBinContent(j, TMath::Sqrt(TMath::Power(h1->GetBinContent(j),2) + TMath::Power(h2->GetBinContent(j),2) + TMath::Power(h3->GetBinContent(j),2)) );
    h1->SetBinError(j,0);
  }
  return;

}

void HistNameHelper::addInQuadrature4(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){
  for(int j = 1; j<h1->GetNbinsX()+1; j++){
    h1->SetBinContent(j, TMath::Sqrt(TMath::Power(h1->GetBinContent(j),2) + TMath::Power(h2->GetBinContent(j),2) + TMath::Power(h3->GetBinContent(j),2) + TMath::Power(h4->GetBinContent(j),2)) );
    h1->SetBinError(j,0);
  }
  return;

}
#endif

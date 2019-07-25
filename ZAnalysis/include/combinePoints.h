#ifndef COMBINEPOINTS
#define COMBINEPOINTS

#include "TMath.h"
#include "TMatrixD.h"
#include <utility>
#include <vector>
#include <iostream>

//see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.55.7731&rep=rep1&type=pdf

class CombinePoints{

public:

  std::vector< double > combine(double e, double m, std::vector< TMatrixD > V);
  TMatrixD getFullCorrMatrix(double dx, double dy);
  TMatrixD getFullUncorrMatrix(double dx, double dy);
private:

};

TMatrixD CombinePoints::getFullCorrMatrix(double dx, double dy){
  TMatrixD temp = TMatrixD(2,2);
  temp[0][0] = dx*dx;
  temp[1][1] = dy*dy;
  temp[0][1] = dx*dy;
  temp[1][0] = dx*dy;
  return temp;
}

TMatrixD CombinePoints::getFullUncorrMatrix(double dx, double dy){
  TMatrixD temp = TMatrixD(2,2);
  temp[0][0] = dx*dx;
  temp[1][1] = dy*dy;
  return temp;
}

//assume the first matrix is a stat error!
std::vector< double > CombinePoints::combine(double e, double m, std::vector< TMatrixD > V){

  //add all the covariance matrices for each uncertainty
  if( V.size() == 0 ) return std::vector<double>();

  TMatrixD C = V.at(0);
  for(unsigned int i = 1; i<V.size(); i++){
    C += V.at(i);
  }
  //C is convariance

  //I is inverse
  TMatrixD I = C;
  I.Invert();
  double norm = -99999;
  norm = I[0][0] + I[1][0] + I[0][1] + I[1][1];

  double w0 = (I[0][0] + I[0][1])/norm;
  double w1 = (I[1][0] + I[1][1])/norm;

  double mean = w0 * e + w1 * m;
  double var2 = w0*w0*C[0][0] + w0*w1*C[0][1] + w1*w0*C[1][0] + w1*w1*C[1][1];

  //do a weighted sum of the sigma^2 of the stat matrix V[0]
  double stat2 = w0 * w0* (V.at(0))[0][0] + w1 * w1 * (V.at(0))[1][1];
  double syst2 = var2 - stat2;

  std::cout << mean << " " << var2 << " " << stat2 << " " << syst2 << " " << w0 << " " << w1 << std::endl;

  std::vector<double> out;
  out.push_back(mean);
  out.push_back(TMath::Sqrt(stat2));
  out.push_back(TMath::Sqrt(syst2));
  out.push_back(w0);
  out.push_back(w1);
  return out;

}

#endif

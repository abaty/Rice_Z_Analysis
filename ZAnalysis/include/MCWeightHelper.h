#ifndef MCWEIGHTHELPER
#define MCWEIGHTHELPER
#include <vector>


class MCWeightHelper{
  public:

    MCWeightHelper();
    ~MCWeightHelper();

    int getSize();
    int getIndx(int i);

  private:

  std::vector< int > weightMap;

};


int MCWeightHelper::getSize(){
  return weightMap.size();
}

int MCWeightHelper::getIndx( int i){
  return weightMap.at(i);
}

MCWeightHelper::MCWeightHelper(){
  weightMap.push_back(1080);//EPPS16nlo
  weightMap.push_back(1177);//nCTEQ15
  weightMap.push_back(0);//free proton
  weightMap.push_back(1);//scale variations
  weightMap.push_back(2);
  weightMap.push_back(3);
  weightMap.push_back(4);
  weightMap.push_back(6);
  weightMap.push_back(8);
  for(int i = 1081; i<= 1176; i++) weightMap.push_back(i);//EPPS16NLO uncertainty copies
}

MCWeightHelper::~MCWeightHelper(){}
#endif

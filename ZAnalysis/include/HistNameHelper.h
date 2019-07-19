#ifndef HISTHELPER
#define HISTHELPER

#include <string>

class HistNameHelper{
  public:

  HistNameHelper();
  ~HistNameHelper();

  std::vector< std::string > name;

};

HistNameHelper::HistNameHelper(){
  name.push_back("massPeak");
  name.push_back("pT");
  name.push_back("y");
  name.push_back("yield");
}

HistNameHelper::~HistNameHelper(){

}

#endif

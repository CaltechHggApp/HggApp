#ifndef dataSetInfo_h
#define dataSetInfo_h

struct dataSetInfo{
  enum dataTypes : signed int {kBackground=-1, kData=0,kSignal=1};
  TString fileName;
  TString label;
  dataTypes type;
  int Ngen;
  bool isList;
  float xsec;
};

#endif

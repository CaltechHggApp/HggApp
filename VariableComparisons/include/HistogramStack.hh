#ifndef HistogramStack_hh
#define HistogramStack_hh

#include <vector>
#include "TString.h"
#include <memory>


template<typename T>
class HistogramStack {
public:
  HistogramStack();
  void Add(T& h);
  T* getTotal(){return hists.back();}
  void Draw(TString options);
private:
  std::vector<T*> hists;

  TString rand_tag;
};



#endif

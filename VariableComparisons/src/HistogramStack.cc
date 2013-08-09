#include "include/HistogramStack.hh"
#include <random>

template<typename T>
HistogramStack<T>::HistogramStack() {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_int_distribution dist();

    rand_tag = TString(dist(re));
}

template<typename T>
void HistogramStack<T>::Add(T& h) {
  T* hclone = h.Clone(h.GetName()+"__HistogramStack_internal_clone_"+rand_tag);
  std::cout << hclone->GetName() << std::endl;
  
  hclone->Add(hists.back());
  hists.push_back(hclone);
}

template<typename T>
void HistogramStack<T>::Draw(TString options) {
  for(auto it = hists.crbegin(); it != hists.crend(); ++it) { //reverse iterate through histograms
    TString this_options = options;
    if(it!=hists.crbegin()) this_options+="SAME";
    it->Draw(this_options);
  }


}

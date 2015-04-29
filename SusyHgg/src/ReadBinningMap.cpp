#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <string>

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <iostream>

std::unique_ptr<std::map<TString,std::vector<int> > > ReadBinningMap(const char* name) {
  ifstream f(name);

  std::unique_ptr<std::map<TString,std::vector<int> > > binmap( new std::map<TString,std::vector<int> > );

  const size_t buffSize=200;
  char buffer[buffSize];
  while( f.getline(buffer,buffSize) ) {
    TObjArray* tok = TString(buffer).Tokenize(" ");
    if(tok->GetEntries()>1) {
      TString key = tok->At(0)->GetName();
      for(int i=1;i<tok->GetEntries();i++) {
	(*binmap)[key].push_back(atoi( tok->At(i)->GetName()));
      }
    }
  }
    
  for(auto& k: *binmap) { 
    std::cout << k.first << " ";
    for(auto v: k.second) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  return binmap;

}

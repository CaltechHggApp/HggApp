#ifndef WeakClassifiers_hh
#define WeakClassifiers_hh

#include "BaseSelector.hh"

#include <string>
#include <vector>
#include <fstream>

class WeakClassifiers : public BaseSelector {
public:
  WeakClassifiers(std::vector<std::string>& fNames, std::string treeName,std::string outputFileName);

private:
  std::fstream outputFile;

  virtual void init();
  virtual void processEntry(Long64_t iEntry);
  virtual void setupOutputTree(){}

};

#endif

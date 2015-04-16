#include "WeakClassifiers.hh"

WeakClassifiers::WeakClassifiers(std::vector<std::string>& fNames, std::string treeName, std::string outputFileName) :
  BaseSelector(),
  outputFile(outputFileName.c_str(),std::fstream::out)
{
  this->loadChain(fNames,treeName);
  this->setDoFill(false);
}


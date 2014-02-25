#include "JECUReader.hh"

#include "fstream"
#include <iostream>

JECUReader::~JECUReader() {
  for(auto it: correctionsArray) {
    delete it.second;
  }
}

bool JECUReader::setCorrections(const char* fileName) {
  std::fstream inFileStream(fileName);
  std::cout << fileName << std::endl;

  if(!inFileStream.is_open()) {
    std::cout << "JECUReader ERROR: Cannot open correction file " << fileName << std::endl;
    return false;
  }
  inFileStream.ignore(1e6,'\n'); //ignore the first line (header)
  
  float e1,e2,e3; //the file is a formatted set of floats.
  //etaMin etaMax 132 ptMin up down ptMin up down ptMin up down ...

  int nRows=0;  //only fill the pt vector on the first row
  
  while( inFileStream >> e1 >> e2 >> e3 ) {


    if(e3==132.) { //e1 and e2 refer to eta points
      correctionsArray.push_back( std::make_pair(e1, new corrVector) ); //build a new entry in the correctionsArray
      nRows++;
    }else{
      if(nRows==1) {
	ptThresholds.push_back(e1);
	//std::cout << e1 << "  " << e2 << std::endl;
      }
      correctionsArray.back().second->push_back( std::make_pair(e2,e3) ); //put the corrections into the vector
    }
  }
  correctionsArray.push_back( std::make_pair(e2,new corrVector) ); //for the last line, push back an empty vector and the upper eta threshold
  //so that we can identify out-of-range jets on the positive side


  int nPtEntries = ptThresholds.size();
  float last_eta=-6;
  //check that everything is the same size
  for(int i=0;i<correctionsArray.size()-1;i++) { //use an index so that we know the line number
    if( correctionsArray.at(i).first < last_eta) {
      std::cout << "JECUReader ERROR: invalid correction (line "<< i << "). Eta values not in ascending order!" << std::endl;
      return false;
    }
    last_eta = correctionsArray.at(i).first;
    if( correctionsArray.at(i).second->size() != nPtEntries ) {
      std::cout << "JECUReader ERROR: invalid correction (line "<<i<<").  number of corrections (" << correctionsArray.at(i).second->size() 
		<<") does not match number of pt bins (" <<  nPtEntries << ")" << std::endl;
      for(int j=0;j<correctionsArray.at(i).second->size();j++) {
	std::cout <<  correctionsArray.at(i).second->at(j).first << std::endl;
      }
      return false;
    }
  }

  __valid=true;
  return true;
}

float JECUReader::getCorrection(float pt, float eta,direction dir) throw(std::runtime_error) {
  if(!__valid) throw std::runtime_error("JECUReader ERROR: attempted to retreive corrections from invalid JECUReader");

  if(eta < correctionsArray.front().first) throw std::runtime_error("JECUReader ERROR: jet eta value out-of-bounds");
  if(pt < ptThresholds.at(0)) throw std::runtime_error("JECUReader ERROR: jet pt value out-of-bounds");

  auto corrArrayIt = correctionsArray.begin();

  //find the first record with a bigger eta than the input eta and step back one record
  for(;corrArrayIt != correctionsArray.end(); corrArrayIt++) {
    if(corrArrayIt->first > eta) {
      corrArrayIt--;
      break;
    }
  }
  //if we hit the end of the loop without finding anything, it means the eta was too big
  if(corrArrayIt == correctionsArray.end()) throw std::runtime_error("JECUReader ERROR: jet eta value out-of-bounds");

  float index=ptThresholds.size()-1;
  for(int i=0;i<ptThresholds.size();i++) { //find the index from the pt
    if(pt< ptThresholds.at(i)) {
      index=i-1; // when we find a bigger pt threshold, set the index to the previous i
      break;
    }
  }//if we never find a bigger pt value, use the corrections for the last index

  if(dir == kUp) {
    return corrArrayIt->second->at(index).first;
  }else{
    return corrArrayIt->second->at(index).second;
  }
}

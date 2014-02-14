/**
   This class will read in a text file from the JetMET POG detailing the JEC uncertainties.
   It then provides a functor that takes in a jet pt and eta and returns the associated uncertainty
*/

#ifndef JECUReader_hh
#define JECUReader_hh

#include <exception>
#include <stdexcept>

#include <vector>
#include <memory>

class JECUReader{
public:
  JECUReader(){} //!< default constructor, correction file must be given elsewhere
  JECUReader(const char* fileName){setCorrections(fileName);} //!< construct and parse the 
  ~JECUReader();

  bool setCorrections(const char* fileName);  //!< sets (or resets) the correction file.  Returns true if successful or false otherwise

  bool isValid(){return __valid;}

  enum direction{kUp,kDown}; //!< defines up or down JEC Uncertainty
  float getCorrection(float pt, float eta, direction dir) throw(std::runtime_error);
private:
  bool __valid=false; //!< set true if the class is initialized correctly.  Trying to run corrections without this set will throw and exception

  /*
    the JEC uncertainties are a matrix of pairs <up,down> where the row is set by the eta of the jet and the column by its pt.
    We will store the map as a vector of <float,const vector*>, where the float is the minimum eta of this of this entry (the maximum is determined
    from the next entry) and the vector points to a vector<pair<float,float>>, where each entry of that is a pair of <up,down> corrections.
    We will store a vector where each entry corresponds to the minimum pt, which will determine the index within the vector of corrections to choose.
  */

  std::vector<float> ptThresholds; //holds the minimum value of pt for this index

  typedef std::vector<std::pair<float,float> > corrVector;

  std::vector< std::pair<float, corrVector* > > correctionsArray;
};

#endif

//-------------------------------------------------------
// Description:
//    Class for SUSY search analyses
// Authors:
//
//-------------------------------------------------------

/// The VertexProjector class takes as input a point in
/// the 3D space and provided an input four momentum
/// vector , measured w.r.t. (0,0,0), calculates the
/// four momentum w.r.t the new vertex. The calculation
/// is done assuming the input object to be massless,
/// i.e. E and p do not change, even if the px, py, and pz
/// do. The class can be used to obatin back the new 
/// coordinates, as well as the missing energy and the
/// missing momentum.

#ifndef VertexProjector_h
#define VertexProjector_h

class VertexProjector {
public:

  /// Class Constructor
  VertexProjector(double x0, double y0, double z0); 
  /// Class Destructor
  virtual ~VertexProjector();     

  /// Set the inner radious of the Ecal
  void SetEcalRad(double rad) {EcalRad = rad;} 

  /// Provide the four momentum to be projected to the new point
  void SetLorentzVector(double thisE, double thispx, double thispy, double thispz);
  
  double GetEt() {return E*sin(Theta);}     ///< Returns the transverse energy
  double GetEtX(){return E*sin(Theta)*cos(Phi);} ///< Returns the x component of the transverse energy
  double GetEtY(){return E*sin(Theta)*sin(Phi);} ///< Returns the y component of the transverse energy
  double GetEta(){return -log(tan(Theta/2.));} ///< Returns the pseudo-rapidity

private:
  /// Calculate the new vector
  void CalcNewVector();
  
  double VtxX0; ///< X coordinate of the new vertex
  double VtxY0; ///< Y coordinate of the new vertex
  double VtxZ0; ///< Z coordinate of the new vertex

  double E;     ///< Energy
  double px;    ///< X component of the momentum
  double py;    ///< Y component of the momentum
  double pz;    ///< Z component of the momentum
  double Theta; ///< Polar angle
  double Phi;   ///< Azimuthal angle
  double EcalRad; ///< The inner radious of Ecal

};
#endif

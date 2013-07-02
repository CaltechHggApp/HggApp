#ifndef _plane_hh_included_
#define _plane_hh_included_
#include <TVector3.h>
#include <TVector2.h>

// Class Plane cooked on the fly. Somebody has better ideas?
class Plane {
public:
  // Default constructor.
  Plane() {}

  // Constructor that takes one vector to define the plane normal.
  Plane(TVector3& n): v1(n.Y(),-n.X(),0) {
    v1 = v1.Unit();
    v2 = n.Cross(v1);
    v2 = v2.Unit();
  }

  // Constructor that takes two vectors to define a plane.
  // These two vectors are orthonormalized (rotating the 2nd one)
  // and are used to define a coordinate system:
  // v1 = (1,0,0), v2 = (0,1,0).
  Plane(TVector3& vec1, TVector3& vec2): v1(vec1), v2(vec2) {
    v1 = v1.Unit();
    v2 = v2 - (v1*v2)*v1;
    v2 = v2.Unit();
  }
  
  // Function to return the normal (defined as v1^v2).
  // In the internal frame of reference, this is n = (0,0,1).
  TVector3 normalVector() {
    return v1.Cross(v2);
  }
  
  // Function to project a vector into the plane. 
  // Returns a 3D vector.
  TVector3 projection(const TVector3& vector) {
    TVector3 theProjectedVector = vector*v1 + vector*v2;
    return theProjectedVector;
  }

  // Function to project a vector into the plane, rotated to the frame
  // where v1 = (1,0,0), v2 = (0,1,0) and the normal is (0,0,1).
  // Returns a 2D vector.
  TVector2 projection2d(const TVector3& vector) {
    TVector2 theProjectedVector (vector*v1,vector*v2);
    return theProjectedVector;
  }

private:
  TVector3 v1;
  TVector3 v2;
  
};
#endif

#include <iostream>
#include "SimpleMatrix.h"

int main(int argc, char** argv){
  SimpleMatrix<double> a(2,2);
  a.set(0,0,1.1);
  a.set(0,1,2.3);
  a.set(1,0,-2.3);
  a.set(1,1,1.2);

  SimpleMatrix<double> b(2,1);
  b.set(0,0,0.5);
  b.set(1,0,-0.5);
  
  SimpleMatrix<double>* c;
  a.dot(c,b);

  std::cout << c->get(0,0) << "   " << c->get(1,0) << std::endl;
  return 0;
}

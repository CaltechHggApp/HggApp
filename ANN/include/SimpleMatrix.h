#ifndef SimpleMatrix_h
#define SimpleMatrix_h

#include <iostream>
#include <cstdlib>
#include <memory>
#include <exception>

template<typename T>
class SimpleMatrix {
public:
  SimpleMatrix(unsigned int r, unsigned int c);
  ~SimpleMatrix();
  inline unsigned int getRows(){ return rows; }
  inline unsigned int getCols(){ return cols; }

  void set(unsigned int i, unsigned int j, T v);
  T get (unsigned int i, unsigned int j);

  void dot(SimpleMatrix *out, SimpleMatrix &right);
private:
  unsigned int rows;
  unsigned int cols;
  T** data=0;
  void clearData();
  bool check(unsigned int i, unsigned int j);
};

#endif

#include <SimpleMatrix.h>
#include <assert.h>
#include <iostream>

template<typename T>
SimpleMatrix<T>::SimpleMatrix(unsigned int r, unsigned int c):
  rows(r),
  cols(c)
{
  if(rows*cols > 100000){
    return;
  }
  try{
    data = new T*[rows];
    for(int i=0;i<rows;i++){
      data[i] = new T[cols];
    }
  }catch(std::exception e){
    std::cout << "Error allocating array of size " << rows << " x " << cols << std::endl;
    std::cout << e.what() << std::endl;
    clearData();
    throw e;
  }
}

template<typename T>
SimpleMatrix<T>::~SimpleMatrix(){
  clearData();
}

template<typename T>
void SimpleMatrix<T>::clearData(){
  for(int i=0;i<rows;i++)
    delete [] data[i];
  delete [] data;
}

template<typename T>
bool SimpleMatrix<T>::check(unsigned int i, unsigned int j){
  assert(i < rows && j < cols);
}

template<typename T>
void SimpleMatrix<T>::set(unsigned int i, unsigned int j, T v){
  check(i,j);
  data[i][j] = v;
}

template<typename T>
T SimpleMatrix<T>::get(unsigned int i, unsigned int j){
  check(i,j);
  return data[i][j];
}

template<typename T>
void SimpleMatrix<T>::dot(SimpleMatrix* out, SimpleMatrix &right){
  out = new SimpleMatrix(rows,right.getCols());

  for(int oR=0;oR<rows;oR++){
    for(int oC=0;oC<right.getCols();oC++){
      T val = 0;
      for(int mC=0;mC<cols;mC++){
	for(int mR=0;mR<right.getRows();mR++){
	  val+=this->get(oR,mC)*right.get(mR,oC);
	}
      }
      out->set(oR,oC,val);
    }
  }
}


//instantiate
template class SimpleMatrix<double>;

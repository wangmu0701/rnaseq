//*********************************************************
// This file is part of Rapsodia released under the LGPL. *
// The full COPYRIGHT notice can be found in the top      *
// level directory of the Rapsodia distribution           *
//*********************************************************
#include <vector>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include "RAinclude.ipp"
#include "HigherOrderTensor.hpp"

extern void readData();
extern void init(const std::vector<RAfloatD>&);
extern void update_eu();
extern RAfloatD likelihood_func();

int main() { 
  int rc=0;
  double myEps=1.0E-12;
  unsigned short  n,o;
  n=7;
  o=2;	
  HigherOrderTensor T(n,o); 
  int dirs=T.getDirectionCount();
  std::cout << "Number of directions: " << dirs << std::endl;
  int i,j,k;
  // argument values
  std::vector<RAfloatD> x(n);
  std::vector<double> xv(n);
  xv[0] = 0.05; // ll0=0.05;
  xv[1] = 7.0;  // l0[0][0]=7.0;
  xv[2] = 8.0;  // l0[1][0]=8.0;
  xv[3] = 9.0;  // l0[2][0]=9.0;
  xv[4] = 0.4;  // p0[0][0]=0.4;
  xv[5] = 0.3;  // p0[1][0]=0.3;
  xv[6] = 0.2;  // p0[2][0]=0.2;

  for(i=1;i<=n;i++) {
    x[i-1] = xv[i-1];
  }

  // get the seed matrix
  Matrix<unsigned int> SeedMatrix=T.getSeedMatrix();
  for(i=1;i<=dirs;i++) { 
    for(j=1;j<=n;j++) { 
      x[j-1].set(i,1,SeedMatrix[j-1][i-1]);
    }
  }
 
  readData();
  // compute the target function
  struct timeval tv1, tv2;
  double time_elapsed = 0;
  gettimeofday(&tv1, NULL);
  init(x);
  update_eu();
  RAfloatD y = likelihood_func();
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  printf("likelihood = %.10f, time elapsed = %.10f\n", y.v, time_elapsed);
  // transfer the taylor coefficients
  Matrix<double> TaylorCoefficients(o,dirs);
  for(i=1;i<=o;i++) { 
    for(j=1;j<=dirs;j++) { 
      TaylorCoefficients[i-1][j-1]=y.get(j,i);
    }
  }
  T.setTaylorCoefficients(TaylorCoefficients);
  // harvest the compressedTensor
  for (k=1; k<=o; k++){ 
    std::cout << "order: " << k << std::endl;
    double entry=0.0;
    std::vector<double> compressedTensor=T.getCompressedTensor(k);
    HigherOrderTensor Helper(n,k); 
    dirs=Helper.getDirectionCount();
    // get the helper seed matrix
    Matrix<unsigned int> HelperSeedMatrix=Helper.getSeedMatrix();
    for(i=1;i<=dirs;i++) {
      std::vector<unsigned short> index(n);
      for(j=1;j<=n;j++) { 
	index[j-1]=HelperSeedMatrix[j-1][i-1];
      }
      std::cout << "T";
      for(j=1;j<=n;j++) { 
        std::cout << "[" << std::setw(2) << HelperSeedMatrix[j-1][i-1] << "]";
      }
      std::cout << " = " << compressedTensor[i-1];
      std::cout << std::endl;
    }
  }
  return rc;
}


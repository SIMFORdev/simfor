#include <iostream>
#include <cmath>
#include <simfor/Vineyard_omp.hpp>

using namespace std;

namespace simfor{

vec Vineyard_omp(vec matr1,int n1, int m1, vec  matr2,  int n2, int m2) {

  vec mulH(n1, 0);
  vec mulV(m2, 0);

  boost::numeric::ublas::vector<double>   res(n1*m2,0);
  #pragma omp parallel for
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < m1 / 2; j++) {
      mulH(i) += matr1(i * m1 + j * 2) * matr1(i * m1 + j * 2 + 1);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < m2; i++) {
    for (int j = 0; j < n2 / 2; j++) {
      mulV(i) += matr2(2*j * m2 + i) * matr2((j * 2 + 1)*m2+i);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < m2; j++) {
      res(i*m2+j) = -mulH(i) - mulV(j);
      for (int k = 0; k < m1 / 2; k++) {
        res(i*m2+j) += (matr1(i * m1 + 2 * k + 1) + matr2(2 * k * m2 + j)) *
                     (matr1(i * m1 + 2 * k) + matr2((2 * k + 1) * m2 + j));
      }
    }
  }

  if (m1 % 2 == 1) {
    #pragma omp parallel for
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < m2; j++) {
        res(i*m2+j) += matr1(i * m1 + m1 - 1) * matr2((m1 - 1) * m2 + j * n2);
      }
    }
  }

  return res;
}
}

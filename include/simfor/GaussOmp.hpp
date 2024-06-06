#ifndef SIMFOR_GAUSS_OMP_HPP_
#define SIMFOR_GAUSS_OMP_HPP_

#include <simfor/types.hpp>

namespace simfor{
    vec GaussianEliminationOmp(matr &mat, int N);
    void SwapRow(matr &mat, int i, int j, int N);
    int ForwardElimOmp(matr &mat, int N);
    vec BackSubOmp(matr &mat, int N);
}

#endif
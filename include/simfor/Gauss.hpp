#ifndef SIMFOR_GAUSS_HPP_
#define SIMFOR_GAUSS_HPP_

#include <simfor/internal/types.hpp>

namespace simfor{
    vec GaussianElimination(matr &mat, int N);
    void SwapRow(matr &mat, int i, int j, int N);
    int ForwardElim(matr &mat, int N);
    vec BackSub(matr &mat, int N);
}

#endif
#ifndef SIMFOR_GRADIENTS_HPP_
#define SIMFOR_GRADIENTS_HPP_

#include <simfor/internal/types.hpp>

namespace simfor{

    vec matrixMulVector(const matr &A, const vec &V);
    vec vectorCombination(double a, const vec &U, double b, const vec &V);
    double innerProduct(const vec &U, const vec &V);
    double vectorNorm(const vec &V);
    //main funcs
    vec ConjugateGradientSolver(const matr &A, const vec &B, size_t iter_num, double eps);
    vec SteepestDescentSolver(const matr &A, const vec &B, size_t iter_num, double eps);
}

#endif
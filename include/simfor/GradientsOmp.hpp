#ifndef SIMFOR_GRADIENTS_OMP_HPP_
#define SIMFOR_GRADIENTS_OMP_HPP_

#include <simfor/types.hpp>

namespace simfor{

    vec matrixMulVector(const matr &A, const vec &V);
    vec vectorCombination(double a, const vec &U, double b, const vec &V);
    double innerProduct(const vec &U, const vec &V);
    double vectorNorm(const vec &V);
    //main funcs
    vec ConjugateGradientSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps);
    vec SteepestDescentSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps);
}

#endif
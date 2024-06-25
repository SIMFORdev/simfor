#ifndef SIMFOR_GRADIENTS_MPI_HPP_
#define SIMFOR_GRADIENTS_MPI_HPP_

#include <simfor/types.hpp>
#include <boost/mpi.hpp>
#include <vector>
namespace mpi = boost::mpi;

namespace simfor{
    vec matrixMulVectorMpi(const matr &A, const vec &V);
    vec vectorCombinationMpi(double a, const vec &U, double b, const vec &V);
    double innerProduct(const vec &U, const vec &V);
    double vectorNorm(const vec &V);
    //main funcs
    vec ConjugateGradientSolverMpi(const matr &A, const vec &B, size_t iter_num, double eps);
    vec SteepestDescentSolverMpi(const matr &A, const vec &B, size_t iter_num, double eps);
}

#endif
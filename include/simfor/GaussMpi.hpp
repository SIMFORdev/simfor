#ifndef SIMFOR_GAUSS_MPI_HPP_
#define SIMFOR_GAUSS_MPI_HPP_

#include <simfor/internal/types.hpp>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

namespace simfor{
    void GaussianEliminationMpi(matr &mat, vec &B, vec &res, int N);
    matr findInvMatGaussJordan(matr &mat_orig, int order);
    void matmult(matr &A, matr &B, matr &C, int n);
}

#endif
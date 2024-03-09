#ifndef SIMFOR_TRIDIAGONAL_MPI_HPP_
#define SIMFOR_TRIDIAGONAL_MPI_HPP_

#include <simfor/internal/types.hpp>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

namespace simfor{
    vec TridiagonalMpi(matr &mat, vec &vec);
}

#endif
#ifndef SIMFOR_ZEIDEL_MPI_HPP_
#define SIMFOR_ZEIDEL_MPI_HPP_

#include <simfor/internal/types.hpp>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;
namespace simfor{
    using vec = boost::numeric::ublas::vector<double>;
    using matr = boost::numeric::ublas::matrix<double>;
    vec ZeidelMpi(matr &mat, vec &vecB, int N);
}

#endif
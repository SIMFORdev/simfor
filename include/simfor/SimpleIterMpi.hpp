#ifndef SIMFOR_SIMPLEITER_MPI_HPP_
#define SIMFOR_SIMPLEITER_MPI_HPP_

#include <simfor/internal/types.hpp>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;
namespace simfor{
    vec SimpleIterMpi(matr &mat, vec &vecB, int N);
}

#endif
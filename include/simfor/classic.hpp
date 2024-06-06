#ifndef SIMFOR_CLASSIC_HPP
#define SIMFOR_CLASSIC_HPP

#include <omp.h>
#include <simfor/types.hpp>

namespace simfor{
    double scalar_mult(vec a, vec b);

    double scalar_mult_omp(vec a, vec b);

    double scalar_mult_mpi(vec a, vec b);
}

#endif

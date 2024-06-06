#ifndef SIMFOR_CLASSIC_HPP_
#define SIMFOR_CLASSIC_HPP_

#include <simfor/types.hpp>

namespace simfor{
	vec classic(vec a, vec b, int n);

    double scalar_mult(vec a, vec b);

    double scalar_mult_omp(vec a, vec b);

    double scalar_mult_mpi(vec a, vec b);
}

#endif

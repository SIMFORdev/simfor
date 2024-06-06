#ifndef SIMFOR_STRASSEN_RECURSIVE_OMP_HPP_
#define SIMFOR_STRASSEN_RECURSIVE_OMP_HPP_

#include <simfor/types.hpp>

#include "simfor/classic_omp.hpp"
#include "simfor/strassen_recursive.hpp"

namespace simfor{
 	vec strassen_recursive_omp(vec a, vec  b, int n);
}

#endif

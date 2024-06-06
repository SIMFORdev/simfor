#ifndef SIMFOR_STRASSEN_OMP_HPP_
#define SIMFOR_STRASSEN_OMP_HPP_

#include <simfor/types.hpp>
#include "simfor/classic_omp.hpp"
#include "simfor/classic.hpp"

namespace simfor{
     vec strassen_omp(vec a,vec b, int n);
}

#endif

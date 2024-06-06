#ifndef SIMFOR_CBLOCK_OMP_HPP_
#define SIMFOR_CBLOCK_OMP_HPP_

#include <simfor/types.hpp>
#include "simfor/classic_omp.hpp"

namespace simfor{
     vec cblock_omp(vec c,vec  b, int n, int N);
}

#endif

#ifndef SIMFOR_CBLOCK_HPP_
#define SIMFOR_CBLOCK_HPP_

#include <simfor/types.hpp>
#include <simfor/classic.hpp>

namespace simfor{
     vec cblock(vec c,vec  b, int n, int N);
     vec cblock_omp(vec c,vec  b, int n, int N);
}

#endif

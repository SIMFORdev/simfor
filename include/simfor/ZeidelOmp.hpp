#ifndef SIMFOR_ZEIDEL_OMP_HPP_
#define SIMFOR_ZEIDEL_OMP_HPP_

#include <simfor/types.hpp>

namespace simfor{
    vec ZeidelOmp(matr &mat, vec &vecB, int N);
}

#endif
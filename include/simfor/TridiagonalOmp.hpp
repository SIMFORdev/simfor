#ifndef SIMFOR_TRIDIAGONAL_OMP_HPP_
#define SIMFOR_TRIDIAGONAL_OMP_HPP_

#include <simfor/internal/types.hpp>

namespace simfor{
    vec TridiagonalOmp(matr &mat, vec &vec);
}

#endif
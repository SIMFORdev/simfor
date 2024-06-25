#ifndef SIMFOR_ELEMENTARY_HPP
#define SIMFOR_ELEMENTARY_HPP

#include <simfor/types.hpp>
#include <omp.h>
#include <boost/mpi.hpp>

namespace simfor {

double v_norm2(vec a, vec b);

double v_norm2_omp(vec a, vec b);

double m_norm4(matr A);

double m_norm4_omp(matr A);

double psum(vec a);

double psum_omp(vec a);

} // namespace simfor

#endif

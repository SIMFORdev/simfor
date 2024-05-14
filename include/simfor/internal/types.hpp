//
// Created by vadim on 19.09.23.
//

#ifndef SIMFOR_TYPES_HPP
#define SIMFOR_TYPES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace simfor {
using matr = boost::numeric::ublas::matrix<double>;
using vec = boost::numeric::ublas::vector<double>;
using scalar_vector = boost::numeric::ublas::scalar_vector<double>;
}

#endif //SIMFOR_TYPES_HPP

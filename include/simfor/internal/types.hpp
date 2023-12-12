//
// Created by vadim on 19.09.23.
//

#ifndef SIMFOR_TYPES_HPP
#define SIMFOR_TYPES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>

namespace simfor {
    using vec = boost::numeric::ublas::vector<double>;
    using matr = boost::numeric::ublas::matrix<double>;
}

#endif //SIMFOR_TYPES_HPP

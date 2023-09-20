//
// Created by vadim on 19.09.23.
//
#include <iostream>
#include "simfor/SomeClass.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

int main(int argc, char **argv) {
	std::cout << "Hello examples in SIMFOR\n";
	simfor::SomeClass testclass;
	testclass.HelloWorld();

	using namespace boost::numeric::ublas;
	matrix<double> m (3, 3);
	for (unsigned i = 0; i < m.size1 (); ++ i)
		for (unsigned j = 0; j < m.size2 (); ++ j)
			m (i, j) = 3 * i + j;
	std::cout << m << std::endl;
	return 0;
}
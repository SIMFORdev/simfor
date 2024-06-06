#include <iostream>
#include <cmath>
#include "simfor/multiplyVectorMatrix.hpp"

using namespace std;

namespace simfor{

vec multiplyVectorMatrix(vec A,vec b,int n){
vec v(n);
for (unsigned i = 0; i < n; i++)
    {
        v(i) = 0;
        for (unsigned j = 0; j < n; j++)
            v(i) += A(i*n+j) * b(j);
    }
    return v;
}
}



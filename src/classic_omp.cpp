 #include <iostream>
#include <cmath>
#include <simfor/classic_omp.hpp>

using namespace std;

namespace simfor{

     vec classic_omp( vec a, vec b, int n){   
     vec c(n*n);
     int i,j,k;
	for(unsigned i = 0;i<n*n;i++)
    c(i)=0;
    #pragma omp parallel for num_threads(4) shared(a,b,c) private(i, j, k)
    for(i = 0;i<n;i++)
       for(k=0;k<n;k++)
           for(j=0;j<n;j++)
              c(i*n+j)+=a(i*n+k)*b(k*n+j);
    return c;
}
}
 

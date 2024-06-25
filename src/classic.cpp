#include <iostream>
#include <cmath>
#include <simfor/classic.hpp>

using namespace std;

namespace simfor{

     vec classic(vec a, vec b, int n){
       vec c(n*n);
	for(unsigned i = 0;i<n*n;i++)
    		c(i)=0;
	for(unsigned i = 0;i<n;i++)
	for(unsigned j=0;j<n;j++)
	for(unsigned k=0;k<n;k++)
		c(k*n+i)+=a(k*n+j)*b(j*n+i);

    return c;
}
}

#include <iostream>
#include <cmath>
#include "simfor/cblock.hpp"
#include "simfor/classic_omp.hpp"

using namespace std;

namespace simfor{

vec cblock_omp(vec c,vec  b, int n, int N)
{
	int i,j,k,ii,jj,kk;
    vec a(n*n,0.0);
    double t=clock();
   if (N>n) a=simfor::classic_omp(c,b,n);
    else{
    #pragma omp parallel for num_threads(4) shared(a,b,c) private(i, j, k,ii,jj,kk)
   for(k=0;k<n;k+=N)
        for(i=0;i<n;i+=N)
            for(j=0;j<n;j+=N)
                for(kk=k;kk<k+N;kk++)
                    for(ii=i;ii<i+N;ii++)
                        for(jj=j;jj<j+N;jj++)
                            a(ii*n+jj)+=c(ii*n+kk)*b(kk*n+jj);

	t=clock()-t;}
	return a;                                  
}
}

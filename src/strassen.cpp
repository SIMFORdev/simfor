#include <iostream>
#include <cmath>
#include "simfor/classic.hpp"
#include "simfor/strassen.hpp"

using namespace std;

namespace simfor{

vec strassen(vec a,vec  b, int n)
{
    vec c(n*n);
    double t;
	if (n<128) c=simfor::classic(a,b,n);
	else
	{
    int m=n;
	n/=2;
	
	 vec  a11(n*n),a12(n*n),a21(n*n),a22(n*n);
	 vec b11(n*n),b12(n*n),b21(n*n),b22(n*n);
	 vec  m1(n*n),m2(n*n),m3(n*n),m4(n*n),m5(n*n),m6(n*n),m7(n*n);
	 vec  c11(n*n),c12(n*n),c21(n*n),c22(n*n);
	
	for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a11(i*n+j) = a(i*m+j);
            a12(i*n+j) = a(i*m+j+n);
            a21(i*n+j) = a((i+n)*m+j);
            a22(i*n+j) = a((i+n)*m+j+n);

            b11(i*n+j) = b(i*m+j);
            b12(i*n+j) = b(i*m+j+n);
            b21(i*n+j) = b((i+n)*m+j);
            b22(i*n+j) = b((i+n)*m+j+n);
        }
    }
	t=clock();
	noalias(m1)= (classic(a11+a22,b11+b22,n));
	noalias(m2)= (classic(a21+a22,b11,n));
	noalias(m3)= (classic(a11,b12-b22,n));
	noalias(m4)= (classic(a22,b21-b11,n));
	noalias(m5)= (classic(a11+a12,b22,n));
	noalias(m6)= (classic(a21-a11,b11+b12,n));
	noalias(m7)= (classic(a12-a22,b21+b22,n));
	c11=m1+m4-m5+m7;
	c12=m3+m5;
	c21=m2+m4;
	c22=m1+m3-m2+m6;
	t=clock()-t;
	for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c(i*m+j) = c11(i*n+j);
            c(i*m+j+n) = c12(i*n+j);
            c((i+n)*m+j) = c21(i*n+j);
            c((i+n)*m+j+n) = c22(i*n+j);
        }
    }
    }
    a.clear (); b.clear ();
    return c; 
    c.clear();
    
}

}

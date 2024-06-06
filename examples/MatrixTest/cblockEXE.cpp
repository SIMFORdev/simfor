#include "simfor/cblock.hpp"
#include <iostream>

using namespace std;

int main(int argc, char **argv){

    int n=256;
     simfor::vec a(n*n),b(n*n),c(n*n);
    for(unsigned i = 0;i<n;i++)
        for(unsigned j=0;j<n;j++){
            a(i*n+j)=1;
            b(i*n+j)=2;
        }
   double t;
   t=clock();
   c=simfor::cblock(a,b,n,128);
   t=clock()-t;
   std::cout << "Answer: " << [&c](){ for (auto &&i : c){std::cout << i << " ";}; return "\n";}();
}

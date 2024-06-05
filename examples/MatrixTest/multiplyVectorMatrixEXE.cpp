#include "simfor/multiplyVectorMatrix.hpp"
#include <iostream>

using namespace std;

int main(int argc, char **argv){

    int n=8;
     simfor::vec vector(n), matrix(n*n), res(n);
    for(unsigned i = 0;i<n;i++)
        for(unsigned j=0;j<n;j++){
            matrix(i*n+j)=1;
            vector(i)=2;
        }
   double t;
   t=clock();
   res=simfor::multiplyVectorMatrix(matrix,vector,n);
   t=clock()-t;
   std::cout << "Answer: " << [&res](){ for (auto &&i : res){std::cout << i << " ";}; return "\n";}();
}

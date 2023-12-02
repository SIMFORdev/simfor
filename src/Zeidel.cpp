#include <iostream>
#include <cmath>
#include "simfor/Zeidel.hpp"

using namespace std;

namespace simfor{

    vec Zeidel(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N);
        long double eps = 1e-9;
        size_t iter{}, iterLim = 100*N;
        long double sum = {};

        do{
            for(auto i = 0; i < N; i++){
                x_n[i] = vecB[i];
                for(auto j = 0; j < i; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n[i] -= mat(i,j) * x_n[j];
                    };
                }

                for(auto j = i; j < N; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n[i] -= mat(i,j) * x[j];
                    };
                }

                x_n[i] /= mat(i,i);
            }

            int flag = 1;
            for(auto i = 0; i < N-1; i++){
                if (fabsf64x(x_n[i] - x[i]) > eps){
                    flag = 0;
                    break;
                }
            }

            for (auto i = 0; i < N; i++){
                x[i] = x_n[i];
            }

            if (flag) break;

        } while (1);
        
        return x;
    }

}
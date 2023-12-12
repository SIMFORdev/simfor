#include <iostream>
#include <cmath>
#include <omp.h>
#include "simfor/ZeidelOmp.hpp"

using namespace std;

namespace simfor{

    vec ZeidelOmp(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N);
        long double eps = 1e-9, sum{}, tmp{};
        size_t iter = 0, iterLim = 100*N;
        int i{}, j{},  pcnt;

        do{
            #pragma omp taskloop default(shared) private(i) //for shared(mat, x_n, x, vecB) private(i, tmp)
            for(i = 0; i < N; i++){
                // x_n[i] = vecB[i];
                tmp = vecB(i);
                #pragma omp parallel for shared(mat, x_n) private(j) reduction(-:tmp)
                for(j = 0; j < i; j++){
                    if (i != j)
                        tmp -= mat(i,j) * x_n(j);
                        // continue;
                }
                x_n(i) = tmp;

                #pragma omp parallel for shared(mat, x) private(j) reduction(-:tmp)
                for(j = i; j < N; j++){
                    if (i != j)
                        tmp -= mat(i,j) * x(j);
                        // continue;
                }

                x_n(i) = tmp /mat(i,i);
            }

            int flag = 1;
            // #pragma omp parallel for shared(x_n, x) private(i)
            for(i = 0; i < N-1; i++){
                if (fabsf64x(x_n(i) - x(i)) > eps){
                    flag = 0;
                    break;
                }
            }

            if (flag) break;

            #pragma omp parallel for shared(x_n, x) private(i)
            for (i = 0; i < N; i++){
                x(i) = x_n(i);
            }

        } while (1);

        return x;
    }

}
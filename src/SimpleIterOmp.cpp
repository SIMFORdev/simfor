#include "omp.h"
#include "simfor/SimpleIterOmp.hpp"

namespace simfor{

    vec SimpleIterOmp(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N); size_t i{}, j{}, iter{};
        long double eps = 1e-9;
        long double sum{}, x1{}, x2{}, x3{};

        #pragma omp parallel for default(none) shared(N, vecB, mat, x) private(i)
        for (i = 0; i < N; i++){
            x[i] = vecB[i] / mat(i,i);
        }

        //start
        do{
            #pragma omp parallel for default(shared) private(i)
            for(i = 0; i < N; i++){
                x_n[i] = vecB[i] / mat(i,i);
                #pragma omp simd// parallel for default(shared) private(j)
                for(j = 0; j < N; j++){
                    if (i != j){
                        x_n[i] -= mat(i,j) / mat(i,i) * x[j];
                    }
                }
            }

            int flag = 1;

            #pragma omp parallel for default(none) shared(N, x, x_n, eps, flag) private(i)
            for(i = 0; i < N-1; i++){
                if (std::abs(x_n[i] - x[i]) > eps){
                    flag = 0;
                    i = N;
                    // break;
                }
            }

            #pragma omp parallel for default(none) shared(N, x, x_n) private(i)
            for (i = 0; i < N; i++){
                x[i] = x_n[i];
            }

            if (flag) break;

        } while (1);
        
        return x;
    }

}
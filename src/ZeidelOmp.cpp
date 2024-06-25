#include <iostream>
#include <cmath>
#include <omp.h>
#include <simfor/ZeidelOmp.hpp>

using namespace std;

namespace simfor{

    /**
     * Solves a system of linear equations A*x = b using the Zeidel's method.
     *
     * @param mat is the coefficient matrix.
     * @param vecB is the right-hand side vector.
     * @param N is the size of the system.
     * @return the solution vector x.
     */
    vec ZeidelOmp(matr &mat, vec &vecB, int N){
        vec x_n(N), x(N); // solution vectors
        long double eps = 1e-9; // precision
        size_t iter = 0, iterLim = 100*N; // iteration limit
        int i{}, j{}, pcnt; // iterators

        do{
            /*
             * Each thread calculates the new value of x_n[i] in parallel.
             * The shared(mat, x_n, x, vecB) clause ensures that each thread has a copy of
             * the shared variables (mat, x_n, x, vecB).
             * The private(i) clause ensures that each thread has its own copy of the private variable i.
             * The reduction(-:tmp) clause ensures that the reduction operation is thread-safe.
             */
            #pragma omp taskloop default(shared) private(i)
            for(i = 0; i < N; i++){
                // x_n[i] = vecB[i];
                long double tmp = vecB(i);
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

            /*
             * Checks if the precision has been reached.
             * The flag is set to 0 if the precision is not reached,
             * otherwise it's set to 1 (precision reached).
             */
            int flag = 1;
            /*
             * Each thread checks the precision in parallel.
             * The shared(x_n, x) clause ensures that each thread has a copy of
             * the shared variables (x_n, x).
             * The private(i) clause ensures that each thread has its own copy of the private variable i.
             */
            #pragma omp parallel for shared(x_n, x) private(i)
            for(i = 0; i < N-1; i++){
                if (fabsf64x(x_n(i) - x(i)) > eps){
                    flag = 0;
                }
            }

            if (flag) break;

            /*
             * Each thread sets the value of x[i] in parallel.
             * The shared(x_n, x) clause ensures that each thread has a copy of
             * the shared variables (x_n, x).
             * The private(i) clause ensures that each thread has its own copy of the private variable i.
             */
            #pragma omp parallel for shared(x_n, x) private(i)
            for (i = 0; i < N; i++){
                x(i) = x_n(i);
            }

        } while (1);

        return x;
    }

}
#include <iostream>
#include <cmath>
#include "simfor/Zeidel.hpp"

using namespace std;

namespace simfor{

    /**
     * @brief Solves a system of linear equations using the Zeidel method
     * 
     * @param mat Coefficient matrix
     * @param vecB Right-hand side vector
     * @param N Number of equations and unknowns
     * @return Solution vector
     */
    vec Zeidel(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N); // current and previous x-vectors
        long double eps = 1e-9; // tolerance
        size_t iter{}, iterLim = 100*N; // limit of iterations
        long double sum = {}; // sum of the differences between x_n and x

        do{
            // Main loop - calculates current x_n and checks for convergence

            for(auto i = 0; i < N; i++){
                // i-th iteration of the outer loop

                x_n(i) = vecB(i); // set current x_n[i] = b[i]

                for(auto j = 0; j < i; j++){
                    // Calculate x_n[i] without the contribution of x_n[j]
                    // where j < i
                    if (i == j)
                        continue;
                    else {
                        x_n(i) -= mat(i,j) * x_n(j);
                    };
                }

                for(auto j = i; j < N; j++){
                    // Calculate x_n[i] without the contribution of x[j]
                    // where j >= i
                    if (i == j)
                        continue;
                    else {
                        x_n(i) -= mat(i,j) * x(j);
                    };
                }

                x_n[i] /= mat(i,i); // divide x_n[i] by the diagonal element a[i,i]
            }

            int flag = 1; // flag = 1 if x_n and x are the same, else flag = 0
            for(auto i = 0; i < N-1; i++){
                // Check for convergence: if |x_n[i] - x[i]| > eps,
                // set flag = 0 and break
                if (fabsf64x(x_n(i) - x(i)) > eps){
                    flag = 0;
                    break;
                }
            }

            for (auto i = 0; i < N; i++){
                // Update previous x-vector x
                x(i) = x_n(i);
            }

            if (flag) break; // if x_n and x are the same, break the loop

        } while (1); // loop until convergence
        
        return x; // return the solution vector
    }


}
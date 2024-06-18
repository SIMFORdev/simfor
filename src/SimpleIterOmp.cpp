#include <omp.h>
#include <simfor/SimpleIterOmp.hpp>

namespace simfor{

    /**
     * @brief Solves system of linear equations A x = b using Jacobi method
     * @details Jacobi method is an iterative method.
     * Algorithm consists of three main steps: initialization, iteration and convergence check.
     * Initialization step is performed in parallel using OpenMP.
     * In each iteration we calculate new approximated solution x_n using Jacobi formula.
     * After that we check if the solution is converged using convergence check step in parallel.
     * If x_n is converged, we break the loop and return the solution.
     * If x_n is not converged, we go back to the initialization step and repeat the process.
     *
     * @param[in] mat square matrix A of the system
     * @param[in] vecB right-hand side vector b of the system
     * @param[in] N size of the system
     *
     * @return x solution of the system A x = b
     */
    vec SimpleIterOmp(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N); size_t i{}, j{}, iter{};
        long double eps = 1e-9;
        long double sum{}, x1{}, x2{}, x3{};

        // Initialization step
        // Divide work between threads using OpenMP
        #pragma omp parallel for default(none) shared(N, vecB, mat, x) private(i)
        for (i = 0; i < N; i++){
            x[i] = vecB(i) / mat(i,i);
        }

        // Jacobi method iterations
        // Loop until convergence
        do{
            // Calculate new approximated solution x_n in parallel
            #pragma omp parallel for default(shared) private(i)
            for(i = 0; i < N; i++){
                x_n(i) = vecB(i) / mat(i,i);
                // Calculate new x_n in parallel using OpenMP simd pragma
                #pragma omp simd
                for(j = 0; j < N; j++){
                    if (i != j){
                        x_n(i) -= mat(i,j) / mat(i,i) * x(j);
                    }
                }
            }

            // Check convergence in parallel
            int flag = 1;
            #pragma omp parallel for default(none) shared(N, x, x_n, eps, flag) private(i)
            for(i = 0; i < N-1; i++){
                if (std::abs(x_n(i) - x(i)) > eps){
                    flag = 0;
                    i = N;
                    // break;
                }
            }

            // If x_n is converged, break the loop
            if (flag) break;

            // If x_n is not converged, repeat the process
            // Update x with x_n in parallel
            #pragma omp parallel for default(none) shared(N, x, x_n) private(i)
            for (i = 0; i < N; i++){
                x(i) = x_n(i);
            }
        } while (1);
        
        return x;
    }


}
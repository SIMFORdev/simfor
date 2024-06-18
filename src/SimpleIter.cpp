
#include <simfor/SimpleIter.hpp>

namespace simfor{

    /**
     * Solves linear system of equations using the simple iterative method.
     * @param mat coefficient matrix of the system
     * @param vecB right-hand side of the system
     * @param N number of equations in the system
     * @return solution vector of the system
     */
    vec SimpleIter(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N); // current and new estimates of the solution
        double eps = 1e-9; // tolerance for the convergence
        size_t iter{}; // number of iterations of the method
        double sum{}; // sum of the absolute difference between the current and new estimates

        // initial guess
        for (auto i = 0; i < N; i++){
            x(i) = vecB(i) / mat(i,i);
        }

        do{
            // flag = 1;
            iter++; // number of iterations
            sum = 0;
            for(auto i = 0; i < N; i++){
                x_n(i) = vecB(i) / mat(i,i);
                for(auto j = 0; j < N; j++){
                    if (i == j) // diagonal element
                        continue;
                    else {
                        x_n(i) -= mat(i,j) / mat(i,i) * x(j); // lower triangular part of the matrix
                    };
                }
                sum += std::abs(x_n(i) - x(i)); // absolute difference between current and new estimates
            }

            int flag = 1;
            for(auto i = 0; i < N-1; i++){ // check the convergence
                if (sum > eps){ // if the difference is larger than the tolerance
                    flag = 0; // the method did not converge
                    break;
                }
            }

            for (auto i = 0; i < N; i++){
                x(i) = x_n(i); // update the solution
            }

            if (flag) break; // if the method converged break the loop

        } while (1);

        return x;
    }



}
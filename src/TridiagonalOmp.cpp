#include <simfor/TridiagonalOmp.hpp>

namespace simfor {
    /**
     * @brief Solves a tridiagonal system Ax = b using the Thomas algorithm
     * 
     * @param A The coefficient matrix of the system (row-major)
     * @param b The right-hand side vector
     * @return The solution vector x
     */
    vec TridiagonalOmp(matr &A, vec &b){
        size_t N(A.size1()); // Number of equations
        vec result(N); // Result vector
        vec U(N), V(N); // Temporary vectors for the Thomas algorithm
        size_t i{0}; // Loop iterator

        // Initialize the first element
        V(0) = A(0, 1) / (-A(0, 0));
        U(0) = -b(0) / (-A(0, 0));

        // Parallel loop for the first N-1 elements
        // (we don't want to parallelize over the last element, as it is special)
        #pragma omp parallel for private(i)
        for (i = 1; i < N-1; i++){
            // Calculate the i-th element of V
            V(i) = A(i, i+1) / (-A(i,i) - A(i,i-1) * V(i-1));

            // Calculate the i-th element of U
            U(i) = (A(i, i-1)*U(i-1) - b(i)) / (-A(i, i) - A(i,i-1)*V(i-1));
        }

        // Initialize the last element of V
        V(N-1) = 0;
        // Calculate the last element of U
        U(N-1) = (A(N-1,N-2) * U(N-2) - b(N-1)) / (-A(N-1, N-1) - A(N-1, N-2) * V(N-2));

        // Store the result of the last element
        result(N-1) = U(N-1);

        // Parallel loop for the remaining N-1 elements
        #pragma omp parallel for private(i)
        for (i = N-1; i >0; i--){
            // Calculate the (i-1)-th element of the result
            result(i-1) = V(i-1) * result(i) + U(i-1);
        }

        return result;
    }
}
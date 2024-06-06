#include <simfor/Tridiagonal.hpp>

namespace simfor {
    /**
     * Solves the tridiagonal system Ax=b
     * @param A The tridiagonal matrix
     * @param b The right-hand side vector
     * @return The solution vector x
     */
    vec Tridiagonal(matr &A, vec &b){
        size_t N(A.size1()); // Number of equations
        vec result(N); // Solution vector
        vec U(N), V(N); // Temporary vectors U and V

        V(0) = A(0, 1) / (-A(0, 0)); // Forward substitution
        U(0) = -b(0) / (-A(0, 0));

        for (size_t i = 1; i < N-1; i++){
            V(i) = A(i, i+1) / (-A(i,i) - A(i,i-1) * V(i-1));

            U(i) = (A(i, i-1)*U(i-1) - b(i)) / (-A(i, i) - A(i,i-1)*V(i-1));
        }
        
        V(N-1) = 0; // Backward substitution
        U(N-1) = (A(N-1,N-2) * U(N-2) - b(N-1)) / (-A(N-1, N-1) - A(N-1, N-2) * V(N-2));

        result(N-1) = U(N-1);
        for (size_t i = N-1; i >0; i--){
            result(i-1) = V(i-1) * result(i) + U(i-1);
        }

        return result;
    }

}
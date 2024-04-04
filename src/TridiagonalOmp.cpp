#include <simfor/TridiagonalOmp.hpp>

namespace simfor {
    vec TridiagonalOmp(matr &A, vec &b){
        size_t N(A.size1());
        vec result(N);
        vec U(N), V(N);
        size_t i{0};
        V(0) = A(0, 1) / (-A(0, 0));
        U(0) = -b(0) / (-A(0, 0));

        #pragma omp parallel for private(i)
        for (i = 1; i < N-1; i++){
            V(i) = A(i, i+1) / (-A(i,i) - A(i,i-1) * V(i-1));

            U(i) = (A(i, i-1)*U(i-1) - b(i)) / (-A(i, i) - A(i,i-1)*V(i-1));
        }
        
        V(N-1) = 0;
        U(N-1) = (A(N-1,N-2) * U(N-2) - b(N-1)) / (-A(N-1, N-1) - A(N-1, N-2) * V(N-2));

        result(N-1) = U(N-1);
        
        #pragma omp parallel for private(i)
        for (i = N-1; i >0; i--){
            result(i-1) = V(i-1) * result(i) + U(i-1);
        }

        return result;
    }
}
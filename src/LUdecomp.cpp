#include "simfor/LUdecomp.hpp"

namespace simfor{

    int LUPDecompose(matr& A, size_t N, double Tol, vec& P){
        int j{}, imax{}; 
        double maxA{}, absA{};

        for (auto i = 0; i <= N; i++)
            P(i) = i; //Unit permutation matrix, P[N] initialized with N

        for (auto i = 0; i < N; i++) {
            maxA = 0.0;
            imax = i;

            for (auto k = i; k < N; k++)
                if ((absA = fabs(A(k,i))) > maxA) { 
                    maxA = absA;
                    imax = k;
                }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i) {
                j =  P(i);
                P(i) = P(imax);
                P(imax) = j;

                //pivoting rows of A
                SwapRow(A, i, j, N);

                //counting pivots starting from N (for determinant)
                P(N)++;
            }

            for (auto j = i + 1; j < N; j++) {
                A(j,i) /= A(i,i);

                for (auto k = i + 1; k < N; k++)
                    A(j,k) -= A(j,i) * A(i,k);
            }
        }

        return 1;  //decomposition done 
    }

    /* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
    * OUTPUT: x - solution vector of A*x=b
    */
    void LUPSolve(matr& A, vec& P, vec& b, size_t N, vec& x) {
        for (auto i = 0; i < N; i++) {
            x(i) = b(P(i));

            for (auto k = 0; k < i; k++)
                x(i) -= A(i,k) * x(k);
        }

        for (int i = N - 1; i >= 0; i--) {
            for (auto k = i + 1; k < N; k++)
                x(i) -= A(i,k) * x(k);

            x(i) /= A(i,i);
        }
    }

    void SwapRow(matr &mat, int i, int j, size_t N){
        for (auto k=0; k<=N; k++){
            double temp = mat(i,k);
            mat(i,k) = mat(j,k);
            mat(j,k) = temp;
        }
    }

    //Those funcs for test purposes
    /* INPUT: A,P filled in LUPDecompose; N - dimension
    * OUTPUT: IA is the inverse of the initial matrix
    */
    void LUPInvert(matr& A, vec& P, int N, matr& IA) {
    
        for (auto j = 0; j < N; j++) {
            for (auto i = 0; i < N; i++) {
                IA(i,j) = P(i) == j ? 1.0 : 0.0;

                for (auto k = 0; k < i; k++)
                    IA(i,j) -= A(i,k) * IA(k,j);
            }

            for (auto i = N - 1; i >= 0; i--) {
                for (auto k = i + 1; k < N; k++)
                    IA(i,j) -= A(i,k) * IA(k,j);

                IA(i,j) /= A(i,i);
            }
        }
    }

    /* INPUT: A,P filled in LUPDecompose; N - dimension. 
    * OUTPUT: Function returns the determinant of the initial matrix
    */
    double LUPDeterminant(matr& A, vec& P, size_t N) {

        auto det = A(0,0);

        for (auto i = 1; i < N; i++)
            det *= A(i,i);

        return int(P(N) - N) % 2 == 0 ? det : -det;
    }

}
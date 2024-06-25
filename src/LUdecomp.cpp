#include <simfor/LUdecomp.hpp>

namespace simfor{

    /**
     * LU decomposition of a matrix A.
     * 
     * INPUT: A - matrix of size NxN, N - size of A
     * OUTPUT: P - permutation matrix, A - LU decomposed matrix
     * RETURN: 1 if decomposition was successful, 0 otherwise
     * 
     * Algorithm: Gaussian elimination with partial pivoting
     */
    int LUPDecompose(matr& A, size_t N, double Tol, vec& P) {
        int j{}, imax{};  //column and row index of pivot
        double maxA{}, absA{}; //maximum element in column and absolute value

        for (auto i = 0; i <= N; i++)
            P(i) = i; //Unit permutation matrix, P[N] initialized with N

        //Gaussian elimination with partial pivoting
        for (auto i = 0; i < N; i++) {
            maxA = 0.0;
            imax = i;

            //Find pivot element in column i
            for (auto k = i; k < N; k++)
                if ((absA = fabs(A(k,i))) > maxA) { 
                    maxA = absA;
                    imax = k;
                }

            //Check for degeneracy
            if (maxA < Tol) return 0; 

            //Swap rows i and imax
            if (imax != i) {
                j =  P(i);
                P(i) = P(imax);
                P(imax) = j;

                //pivoting rows of A
                SwapRow(A, i, j, N);

                //counting pivots starting from N (for determinant)
                P(N)++;
            }

            //Elimination step
            for (auto j = i + 1; j < N; j++) {
                A(j,i) /= A(i,i);

                for (auto k = i + 1; k < N; k++)
                    A(j,k) -= A(j,i) * A(i,k);
            }
        }

        return 1;  //decomposition done 
    }

    /**
     * Solves the linear system Ax=b using the LUP decomposition of A
     * 
     * INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
     * OUTPUT: x - solution vector of A*x=b
     */
    void LUPSolve(matr& A, vec& P, vec& b, size_t N, vec& x) {
        // Forward substitution step:
        for (auto i = 0; i < N; i++) {
            x(i) = b(P(i)); // Use the initial order in P to start

            // Subtract the effects of all previous rows
            for (auto k = 0; k < i; k++)
                x(i) -= A(i,k) * x(k);
        }

        // Backward substitution step:
        for (int i = N - 1; i >= 0; i--) {
            // Subtract the effects of all later rows
            for (auto k = i + 1; k < N; k++)
                x(i) -= A(i,k) * x(k);

            x(i) /= A(i,i); // Divide by diagonal element to get the solution.
        }
    }

    void SwapRow(matr &mat, int i, int j, size_t N){
        std::swap_ranges(mat.begin1()+i*N, mat.begin1()+(i+1)*N, mat.begin1()+j*N);
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
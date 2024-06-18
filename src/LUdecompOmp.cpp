#include <simfor/LUdecompOmp.hpp>

namespace simfor{

    /**
     * @brief LUP decomposition of a matrix in parallel with OpenMP
     * @details Performs an LUP decomposition of a given matrix A with size N x N,
     *          using partial pivoting. The decomposition is done in place,
     *          so A will be modified. The permutation matrix P is also returned.
     *          The function returns 1 if the decomposition is successful,
     *          or 0 if the matrix is degenerate.
     *
     * @param A Matrix to be decomposed
     * @param N Size of the matrix
     * @param Tol Tolerance below which pivoting is skipped
     * @param P Permutation matrix
     *
     * @return 1 if successful, 0 otherwise
     */
    int LUPDecomposeOmp(matr& A, size_t N, double Tol, vec& P){
        int i{},j{}, k{}, imax{};  //loop counters
        double maxA{}, absA{};     //maximum value in column i, absolute value of A(k,i)

        #pragma omp parallel for private(i)
        for (i = 0; i <= N; i++) //initialize permutation matrix and set P[N] to N
            P(i) = i;

        // #pragma omp parallel for default(shared) private(i, maxA, imax, k, j)
        for (i = 0; i < N; i++) {
            maxA = 0.0;
            imax = i;

            #pragma omp parallel for default(shared) private(k)
            for (k = i; k < N; k++) {
                #pragma omp critical
                if ((absA = fabs(A(k,i))) > maxA) {  //find maximum in column i
                    maxA = absA;
                    imax = k;
                }
            }

            if (maxA < Tol) { //failure, matrix is degenerate
                k = 0;
            }

            if (imax != i) {
                #pragma omp parallel
                #pragma omp taskgroup
                {
                j = P(i);
                P(i) = P(imax);
                P(imax) = j;

                //pivoting rows of A
                SwapRowOmp(A, i, j, N);

                //counting pivots starting from N (for determinant)
                P(N)++;
                }
            }

            #pragma omp parallel for default(shared) private(j, k)
            for (j = i + 1; j < N; j++) {
                A(j,i) /= A(i,i);
                #pragma omp parallel for default(shared) private(k) firstprivate(j)
                for (k = i + 1; k < N; k++)
                    A(j,k) -= A(j,i) * A(i,k);
            }
        }

        return 1;  //decomposition done
    }


    /**
     * @brief Solves linear system A*x=b using LUP decomposition with OpenMP
     * @details Solves the linear system A*x=b where A is a square matrix of
     *          size NxN and b is a vector of length N. The function assumes
     *          that A and P are already filled with the LUP decomposition of
     *          the original matrix A, as computed by LUPDecomposeOmp.
     *          The solution x is computed and stored in the vector x.
     *
     * @param A Matrix A of the linear system
     * @param P Permutation matrix of LUP decomposition
     * @param b Vector b of the linear system
     * @param N Size of the matrix and vectors
     * @param x Solution vector x
     */
    void LUPSolveOmp(matr& A, vec& P, vec& b, size_t N, vec& x) {
        int i{}, k{};  //loop counters

        //forward substitution
        #pragma omp taskloop default(shared) private(i, k)
        for (i = 0; i < N; i++) {
            x(i) = b(P(i));  //use values in original order
            #pragma omp taskloop default(shared) private(k)
            for (k = 0; k < i; k++)
                x[i] -= A(i,k) * x(k);
        }

        //backward substitution
        #pragma omp taskloop default(shared) private(i, k)
        for (i = N - 1; i >= 0; i--) {
            #pragma omp taskloop default(shared) private(k)
            for (k = i + 1; k < N; k++)
                x(i) -= A(i,k) * x(k);

            x(i) /= A(i,i);
        }
    }

    void SwapRowOmp(matr &mat, int i, int j, size_t N){
        size_t k{};
        #pragma omp parallel for private(k) order(concurrent)
        for (k=0; k<=N; k++){
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
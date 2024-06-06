#ifndef SIMFOR_LUDECOMP_HPP_
#define SIMFOR_LUDECOMP_HPP_

#include <simfor/types.hpp>

namespace simfor{
    /* INPUT: A - array of pointers to rows of a square matrix having dimension N
    *        Tol - small tolerance number to detect failure when the matrix is near degenerate
    * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
    *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
    *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
    *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
    */
    int LUPDecompose(matr& A, size_t N, double Tol, vec& P);

    /* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
    * OUTPUT: x - solution vector of A*x=b
    */
    void LUPSolve(matr& A, vec& P, vec& b, size_t N, vec& x);

    /* INPUT: A,P filled in LUPDecompose; N - dimension
    * OUTPUT: IA is the inverse of the initial matrix
    */
    void LUPInvert(matr& A, vec& P, int N, matr& IA);

    /* INPUT: A,P filled in LUPDecompose; N - dimension. 
    * OUTPUT: Function returns the determinant of the initial matrix
    */
    double LUPDeterminant(matr& A, vec& P, size_t N);

    void SwapRow(matr& mat, int i, int j, size_t N);

}
#endif
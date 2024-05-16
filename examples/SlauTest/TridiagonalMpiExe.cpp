#include "simfor/TridiagonalMpi.hpp"

/**
 * @brief Generate a random tridiagonal matrix
 * 
 * The resulting matrix will have the following structure:
 * 
 *       a_ii       a_ij      | i,j \in [0, n-1]
 * a_ij = c_ij 0   0   0    | i == j
 *       c_ij a_i,j 0   0    | i == j+1
 *       0   c_ij a_i,j 0    | i == j-1
 *       0   0   c_ij a_i,j  | otherwise
 * 
 * @param n Number of rows and columns of the generated matrix
 * @return A random tridiagonal matrix with the specified size
 */
simfor::matr genMatNNB(int n){
        simfor::matr m(n, n); // initialize empty matrix
        // Loop through all the elements of the matrix
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                // If the element is on the diagonal
                if (i==j)
                {
                    // Set the element to a random integer between 11 and 110 with
                    // absolute value
                    m(i,j) = 10*fabsf64x(rand()%100+11);
                }
                // If the element is on the sub/super-diagonal
                else if (i == j+1 || i == j-1)
                {
                    // Set the element to a random integer between 0 and 9
                    m(i,j) = rand()%10;
                }
                // Otherwise (the element is not on the diagonal, subdiagonal or
                // superdiagonal)
                else
                {
                    // Set the element to 0
                    m(i,j) = 0;
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector of size n with random values
 * 
 * @param n Size of the vector to generate
 * @return A vector of size n with random values
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); // Initialize a vector of size n
    for(auto i = 0; i < n; /** Loop through all the elements of the vector */
        v[i++] = 10*fabsf64x(rand()%10+11)); /** Set the current element to a random value between 11 and 110 with absolute value */
    return v; /** Return the generated vector */
}

int main(int argc, char *argv[]){
    mpi::environment env (argc, argv);
    mpi::communicator world;

    const int N = 4;
    // simfor::matr myMatA = genMatNNB(N);
    // simfor::vec myVecA = genVecN(N);

    simfor::matr myMatA(N,N);
    simfor::vec myVecA(N);

    //Answer: 1.11859 1.31062 1.50319 1.70798 
    std::vector<std::vector<double>> myMatB = {{ 10.8000, 0.0475,      0, 0     },
                                            {  0.0321, 9.9000, 0.0523, 0     },
                                            {       0, 0.0369, 9.0000, 0.0570},
                                            {       0,      0, 0.0416, 8.1000}};
    std::vector<double> myVecB = {12.1430, 13.0897, 13.6744, 13.8972};

    for (auto i = 0; i < N; i++){
        for (auto j = 0; j < N; j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myVecA(i) = myVecB[i];
    }

    simfor::vec resVec = simfor::TridiagonalMpi(myMatA, myVecA);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}

#include "simfor/Tridiagonal.hpp"

/**
 * @brief Generates a n x n matrix with random values
 *
 * @param n size of the matrix to generate
 *
 * @return simfor::matr n x n matrix with random values
 */
simfor::matr genMatNNB(int n)
{
    simfor::matr m(n, n); // Initialize a n x n matrix

    // Fill the matrix with random values
    // The diagonal elements are between [10, 110]
    // The off-diagonal elements are between [0, 10]
    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            if (i == j) // Diagonal elements
            {
                m(i, j) = 10.0 * fabsf64x(rand() % 100 + 11);
            }
            else if (i == (j + 1)) // Off-diagonal elements to the right
            {
                m(i, j) = rand() % 10;
            }
            else if (i == (j - 1)) // Off-diagonal elements to the left
            {
                m(i, j) = rand() % 10;
            }
            else // Off-diagonal elements to the upper and lower
            {
                m(i, j) = 0;
            }
        }
    }
    return m;
}

/**
 * @brief Generates a vector of size n with random values in [10, 20]
 *
 * @param n size of the vector to generate
 *
 * @return simfor::vec n sized vector with random values in [10, 20]
 */
simfor::vec genVecN(int n)
{
    simfor::vec v(n); // Initialize a n sized vector

    // Fill the vector with random values between 10 and 20
    for (auto i = 0; i < n; v[i++] = 10.0 * fabsf64x(rand() % 10 + 11)); 

    return v;
}


int main(int argc, char const *argv[]){
    
    const int N = 4;
    // simfor::matr myMatA = genMatNNB(N);
    // simfor::vec myVecb = genVecN(N);

    simfor::matr myMatA(N,N);
    simfor::vec myVecb(N);

    //Answer: 1.11859 1.31062 1.50319 1.70798 
    std::vector<std::vector<double>> myMatB = {{ 10.8000, 0.0475,      0, 0     },
                                        {  0.0321, 9.9000, 0.0523, 0     },
                                        {       0, 0.0369, 9.0000, 0.0570},
                                        {       0,      0, 0.0416, 8.1000}};
    std::vector<double> myVecB = {12.1430, 13.0897, 13.6744, 13.8972};

    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myVecb(i) = myVecB[i];
    }

    simfor::vec resVec = simfor::Tridiagonal(myMatA, myVecb);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}

#include "simfor/Zeidel.hpp"

/**
 * @brief Generate matrix of size n x n+1.
 *
 * Matrix is filled with random values, with main diagonal elements
 * having double the value of other diagonal elements.
 *
 * @param n size of matrix
 * @return matrix of size n x n+1
 */
simfor::matr genMatNMB(int n){
        // Matrix of size n x n+1
        simfor::matr m(n, n+1);
        // Iterator for rows
        for(auto i=0;i<n;i++){
            // Iterator for columns
            for(auto j=0;j<n+1;j++){
                // If we are on the main diagonal, set value
                if (i==j)
                {
                    // Value is 10 times the size of the matrix times a random number between 11 and 110
                    m(i,j) = 10*n*fabsf64x(rand()%100+11);
                }
                else
                {
                    // Otherwise set value to a random number between 0 and 9
                    m(i,j) = rand()%10;
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector of size n with random values in [10, 110]
 *
 * @param n size of vector
 * @return simfor::vec vector of size n with random values in [10, 110]
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); // Create a vector of size n
    for(auto i = 0; i < n; /* increment i */) /* Assign value to v[i] and increment i */
        v[i++] = 10*fabsf64x(rand()%10+11); // Set element to 10 * abs(random number) + 11
    return v; // Return the generated vector
}


int main(int argc, char const *argv[])
{
    const auto N = 5;

    //Наша исходная матрица и вектор свободных чисел 
    std::vector<std::vector<double>> myMatB = {{32, 2, 1, 3, 1},
                    {1, 8, 3, 1, 3},
                    {1, 2, 16, 3, 1},   
                    {1, 2, 3, 56, 1},  
                    {1, 2, 1, 3, 32}};
    std::vector<double> myVecB = {43, 14, -3, 169, -13};
    
    //Матрица вида A необходима для работы функции GaussianElimination
    simfor::matr myMatA(N, N);
    simfor::vec myVecA(N);

    // //Копируем значения
    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myVecA(i) = myVecB[i];
    }

    // simfor::matr myMatA = genMatNMB(N);
    // simfor::vec myVecA = genVecN(N);

    simfor::vec resVec = simfor::Zeidel(myMatA, myVecA, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}


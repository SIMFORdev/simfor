#include "simfor/ZeidelOmp.hpp"

/**
 * @brief Generate matrix with main diagonal filled with
 * random numbers and off-diagonal elements filled with
 * random numbers from 0 to 9
 *
 * @param n Size of the matrix
 *
 * @return n x (n+1) matrix
 */
simfor::matr genMatNMB(int n){
        // Matrix of size n x (n+1)
        simfor::matr m(n, n+1);
        // Fill elements
        for(auto i=0;i<n;i++){
            for(auto j=0;j<n+1;j++){
                // If diagonal element
                if (i==j)
                {
                    // Fill with random value from 10 to 100
                    m(i,j) = 10*n*fabsf64x(rand()%100+11);
                }else{
                    // Fill with random value from 0 to 9
                    m(i,j) = rand()%10;
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector of size n with random values
 *
 * @param n size of the vector to generate
 *
 * @return simfor::vec of size n with random values
 */
simfor::vec genVecN(int n)
{
    simfor::vec v(n); // Initialize a vector of size n

    // Fill the vector with random values
    // Each element is a random number between [10, 100]
    for(auto i = 0; i < n; i++)
    {
        v[i] = 10*fabsf64x(rand()%100+11);
    }

    return v;
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

    //Копируем значения
    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myVecA(i) = myVecB[i];
    }

    // simfor::matr myMatA = genMatNMB(N);
    // simfor::vec myVecA = genVecN(N);

    simfor::vec resVec = simfor::ZeidelOmp(myMatA, myVecA, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}
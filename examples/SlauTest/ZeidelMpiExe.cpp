#include "simfor/ZeidelMpi.hpp"

/**
 * @brief Generate matrix of size n x n+1 filled with random numbers
 * 
 * @param n size of matrix
 * @return generated matrix
 */
simfor::matr genMatNMB(int n){
    // matrix of size n x n+1
    simfor::matr m(n, n+1);
    // iterate through rows
    for(auto i=0;i<n;i++){
        // iterate through columns
        for(auto j=0;j<n+1;j++){
            // if diagonal element fill it with random number
            if (i==j)
            {
                // multiply random number by 10*n to make diagonal elements big
                m(i,j) = 10*n*fabsf64x(rand()%100+11);
            }
            // else fill with random number from 0 to 9
            else
            {
                m(i,j) = rand()%10;
            }
        }
    }
    return m;
}

/**
 * @brief Generate a vector of size n filled with random numbers
 * 
 * @param n size of vector
 * @return generated vector
 */
simfor::vec genVecN(int n)
{
    simfor::vec v(n); // Initialize a vector of size n
    for(auto i = 0; i < n; /* increment is done in the body */) 
    {
        // Assign a random number with absolute value between 11 and 20 to the current element
        v[i++] = 10*fabsf64x(rand()%10+11);
    }
    return v;
}

int main(int argc, char *argv[])
{
    const auto N = 5;
    mpi::environment env (argc, argv);
    mpi::communicator world;
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

    // Копируем значения
    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myVecA(i) = myVecB[i];
    }

    // simfor::matr myMatA = genMatNMB(N);
    // simfor::vec myVecA = genVecN(N);

    simfor::vec resVec = simfor::ZeidelMpi(myMatA, myVecA, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}


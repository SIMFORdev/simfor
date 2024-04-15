#include "simfor/GaussOmp.hpp"

/**
 * @brief Generate matrix with random numbers
 *
 * Function generates matrix with N rows and N+1 columns
 * with random numbers. Diagonal elements are integers from
 * 10 to 1000 with absolute value, other elements are integers
 * from 0 to 9.
 *
 * @param n Number of rows and columns of the matrix
 * @return Simfor matrix with random numbers
 */
simfor::matr genMatNMB(int n){
        /* Generate matrix with random numbers */
        simfor::matr m(n, n+1);
        /* Iterate through the matrix */
        for(auto i=0;i<n;i++){
            for(auto j=0;j<n+1;j++){
                /* Set diagonal elements to 10*abs(rand()%100+11) */
                if (i==j)
                {
                    m(i,j) = 10*fabsf64x(rand()%100+11);
                }
                /* Set other elements to rand()%10 */
                else
                {
                    m(i,j) = rand()%10;
                }
            }
        }
        /* Return generated matrix */
        return m;
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
    
    //Матрица вида Ab необходима для работы функции GaussianElimination
    simfor::matr myMatA(N, N+1);

    //Копируем значения
    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMatA(i, j) = myMatB[i][j];
        }
        myMatA(i, N) = myVecB[i];
    }

    // simfor::matr myMatA = genMatNMB(N);

    simfor::vec resVec = simfor::GaussianEliminationOmp(myMatA, N);

    std::cout << "Answer: " << [resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}


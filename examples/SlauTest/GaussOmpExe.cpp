#include "simfor/GaussOmp.hpp"

simfor::matr genMatNMB(int n){
        simfor::matr m(n, n+1);
        for(auto i=0;i<n;i++){
            for(auto j=0;j<n+1;j++){
                if (i==j)
                {
                    m(i,j) = 10*fabsf64x(rand()%100+11);
                }else{
                    m(i,j) = rand()%10;
                }
            }
        }
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


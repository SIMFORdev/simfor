#include "simfor/ZeidelMpi.hpp"

simfor::matr genMatNMB(int n){
        simfor::matr m(n, n+1);
        for(auto i=0;i<n;i++){
            for(auto j=0;j<n+1;j++){
                if (i==j)
                {
                    m(i,j) = 10*n*fabsf64x(rand()%100+11);
                }else{
                    m(i,j) = rand()%10;
                }
            }
        }
        return m;
}

simfor::vec genVecN(int n){
    simfor::vec v(n);
    for(auto i = 0; i < n; v[i++] = 10*fabsf64x(rand()%10+11)); 
    return v;
}

int main(int argc, char *argv[])
{
    const auto N = 8192;
    namespace mt  = mpi::threading;
    mpi::environment env (argc, argv, mt::multiple, 1);
    mpi::communicator world;
    //Наша исходная матрица и вектор свободных чисел 
    std::vector<std::vector<double>> myMatB = {{32, 2, 1, 3, 1},
                    {1, 8, 3, 1, 3},
                    {1, 2, 16, 3, 1},   
                    {1, 2, 3, 56, 1},  
                    {1, 2, 1, 3, 32}};
    std::vector<double> myVecB = {43, 14, -3, 169, -13};
    
    //Матрица вида A необходима для работы функции GaussianElimination
    // simfor::matr myMatA(N, N);
    // simfor::vec myVecA(N);

    //Копируем значения
    // for (auto i = 0; i < myMatB.size(); i++){
    //     for (auto j = 0; j < myMatB[i].size(); j++){
    //         myMatA(i, j) = myMatB[i][j];
    //     }
    //     myVecA(i) = myVecB[i];
    // }

    simfor::matr myMatA = genMatNMB(N);
    simfor::vec myVecA = genVecN(N);

    simfor::vec resVec = simfor::ZeidelMpi(myMatA, myVecA, N);

    // std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}


#include "simfor/TridiagonalMpi.hpp"

simfor::matr genMatNNB(int n){
        simfor::matr m(n, n);
        for(auto i=0;i<n;i++){
            for(auto j=0;j<n;j++){
                if (i==j)
                {
                    m(i,j) = 10*fabsf64x(rand()%100+11);
                }else if (i == (j+1))
                {
                    m(i,j) = rand()%10;
                }else if (i == (j-1))
                {
                    m(i,j) = rand()%10;
                }else{
                    m(i,j) = 0;
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

int main(int argc, char *argv[]){

    namespace mt  = mpi::threading;
    mpi::environment env (argc, argv, mt::multiple, 1);
    mpi::communicator world;

    const int N = 2048;
    matr myMatA = genMatNNB(N);
    vec myVecb = genVecN(N);

    // simfor::matr myMatA(N,N);
    // simfor::vec myVecb(N);

    //Answer: 1.11859 1.31062 1.50319 1.70798 
    std::vector<std::vector<double>> myMatB = {{ 10.8000, 0.0475,      0, 0     },
                                            {  0.0321, 9.9000, 0.0523, 0     },
                                            {       0, 0.0369, 9.0000, 0.0570},
                                            {       0,      0, 0.0416, 8.1000}};
    std::vector<double> myVecB = {12.1430, 13.0897, 13.6744, 13.8972};

    // for (auto i = 0; i < N; i++){
        // for (auto j = 0; j < N; j++){
            // myMatA(i, j) = myMatB[i][j];
        // }
        // myVecb(i) = myVecB[i];
    // }

    simfor::vec resVec = simfor::TridiagonalMpi(myMatA, myVeca);

    // std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}

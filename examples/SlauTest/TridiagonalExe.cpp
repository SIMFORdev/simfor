#include "simfor/Tridiagonal.hpp"

int main(int argc, char const *argv[]){
    
    const int N = 4;

    simfor::matr myMatA(N,N);
    simfor::vec myVecb(N);

    //Answer: 1.11859 1.31062 1.50319 1.70798 
    std::vector<std::vector<double>> myMatB = {{ 10.8000, 0.0475,      0, 0     },
                                        {  0.0321, 9.9000, 0.0523, 0     },
                                        {       0, 0.0369, 9.0000, 0.0570},
                                        {       0,      0, 0.0416, 8.1000}};
    std::vector myVecB = {12.1430, 13.0897, 13.6744, 13.8972};

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

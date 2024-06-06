#include "simfor/TridiagonalOmp.hpp"

/**
 * @brief Generate a random matrix with non-zero elements in the main diagonal and two neighboring diagonals.
 * @param n The number of rows and columns of the matrix.
 * @return A random matrix with non-zero elements in the main diagonal and two neighboring diagonals.
 */
simfor::matr genMatNNB(int n){
        simfor::matr m(n, n); // create a matrix with n rows and columns
        for(auto i=0;i<n;i++){ // iterate over rows
            for(auto j=0;j<n;j++){ // iterate over columns
                if (i==j) // if i == j, it is in the main diagonal
                {
                    m(i,j) = 10*fabsf64x(rand()%100+11); // set the element to a random number in [11,110]
                }else if (i == (j+1)) // if i == j+1, it is in the upper diagonal
                {
                    m(i,j) = rand()%10; // set the element to a random number in [0,9]
                }else if (i == (j-1)) // if i == j-1, it is in the lower diagonal
                {
                    m(i,j) = rand()%10; // set the element to a random number in [0,9]
                }else{
                    m(i,j) = 0; // otherwise the element is 0
                }
            }
        }
        return m;
}


/**
 * @brief Generate a vector of length n with random values in [11,110]
 * @param n The length of the vector to generate
 * @return A vector of length n with random values in [11,110]
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); // initialize vector of length n
    for(auto i = 0; i < n; /* increment index inside loop */) {
        v[i++] = 10*fabsf64x(rand()%10+11); // set value and increment index
    }
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

    simfor::vec resVec = simfor::TridiagonalOmp(myMatA, myVecb);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

    return 0;
}

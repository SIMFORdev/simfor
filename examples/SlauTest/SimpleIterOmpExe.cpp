#include "simfor/SimpleIterOmp.hpp"

/**
 * @brief Generate a matrix with random numbers
 * 
 * @param n Dimension of the matrix
 * @return simfor::matr Generated matrix
 */
simfor::matr genMatNNB(int n) {
    // Initialize matrix with dimension n x n
    simfor::matr m(n, n);

    // Iterate through rows
    for(auto i=0;i<n;i++) {
        // Iterate through columns
        for(auto j=0;j<n;j++) {
            // If diagonal element
            if (i==j) {
                // Assign a random number with absolute value between 11 and 20
                m(i,j) = fabsf64x(rand()%10+11);
            } else {
                // Assign a random number between 0 and 9
                m(i,j) = rand()%10;
            }
        }
    }

    return m;
}

/**
 * @brief Generate a vector with n elements
 * 
 * @param n Size of the vector
 * @return simfor::vec Vector with n elements
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); /* Initialize a vector with n elements */

    /* Iterate over vector and assign random values to each element */
    for(auto i = 0; i < n; /* Increment i */ ) 
    {
        /* Assign a random value to the current element */
        v[i++] = 10*fabsf64x(rand()%10+11); /* [10, 20] */
    }
    return v; /* Return the generated vector */
}

int main(int argc, char** argv){
	const auto N = 5;
    
    // simfor::matr Mat = genMatNNB(N);
    // simfor::matr copyMat = Mat;
    // simfor::matr myMat(N, N);
    // simfor::vec myVec = genVecN(N);

    // for(auto i = 0; i < N; i++)
	// {
	// 	for(auto j = 0; j < N; j++)
	// 	{
	// 		myMat(i,j) =  Mat(j,i) * copyMat(i,j);
	// 	}
	// }

    // Uncomment this code to test the example; The answer is ¬=(1,1.92,-1,3,-0.8)
    // Of course do not forget to comment out code above
    simfor::matr myMat(N, N);
    simfor::vec myVec(N);

    std::vector<std::vector<double>> myMatB = {{32, 2, 1, 3, 1},
                    {1, 8, 3, 1, 3},
                    {1, 2, 16, 3, 1},   
                    {1, 2, 3, 56, 1},  
                    {1, 2, 1, 3, 32}};
    std::vector<double> myVecB = {43, 14, -3, 169, -13};
    for (auto i = 0; i < myMatB.size(); i++){
        for (auto j = 0; j < myMatB[i].size(); j++){
            myMat(i, j) = myMatB[i][j];
        }
        myVec(i) = myVecB[i];
    }

    simfor::vec resVec = simfor::SimpleIterOmp(myMat, myVec, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

	return 0;
}
#include "simfor/SimpleIter.hpp"

/**
 * @brief Generates a square matrix of size n with random values
 * 
 * @param n size of the matrix to generate
 * @return simfor::matr n x n sized matrix with random values
 */
simfor::matr genMatNNB(int n)
{
    simfor::matr m(n, n); // Initialize a n x n matrix

    /* Fill the diagonal with absolute values of random integers
     * between 11 and 20
     */
    for(auto i=0; i<n; i++)
    {
        m(i, i) = fabsf64x(rand() % 10 + 11);
    }

    /* Fill the rest of the matrix with random integers between 0 and 9
     */
    for(auto i=0; i<n; i++)
    {
        for(auto j=0; j<n; j++)
        {
            if (i != j)
            {
                m(i, j) = rand() % 10;
            }
        }
    }

    return m;
}

/**
 * @brief Generates a vector of size n with random values
 * 
 * @param n size of the vector to generate
 * @return simfor::vec n sized vector with random values
 */
simfor::vec genVecN(int n)
{
    simfor::vec v(n); // Initialize a vector of size n
    for (auto i = 0; i < n; /* increment in the for loop */)
    {
        /* Set the i-th value of the vector to a random value
         * between 10 and 20 multiplied by the absolute value
         * of a random integer between 0 and 10
         */
        v[i++] = 10 * fabsf64x(rand() % 10 + 11);
    }
    return v;
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

    simfor::vec resVec = simfor::SimpleIter(myMat, myVec, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

	return 0;
}
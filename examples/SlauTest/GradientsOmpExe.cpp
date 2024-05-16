#include "simfor/GradientsOmp.hpp"

/**
 * @brief Generate matrix with main diagonal equal to 10*|random integer between 11 and 100|
 *        and off-diagonal elements equal to random integers between 0 and 9
 * @param n - size of matrix
 * @return matrix with given properties
 */
simfor::matr genMatNNB(int n)
{
        simfor::matr m(n, n); // create matrix with n columns and n rows
        for (auto i = 0; i < n; i++) // iterate over rows
        {
            for (auto j = 0; j < n; j++) // iterate over columns
            {
                if (i == j) // if current element is on the main diagonal
                {
                    m(i, j) = 10 * fabsf64x(rand() % 100 + 11); // set it to 10 times absolute value of random integer between 11 and 100
                }
                else // if current element is not on the main diagonal
                {
                    m(i, j) = rand() % 10; // set it to random integer between 0 and 9
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector of length n with random values in [10, 100]
 * @param n Length of the generated vector
 * @return Vector with length n filled with random values in [10, 100]
 */
simfor::vec genVecN(int n){
    simfor::vec v(n);  // Create a vector of length n
    for(auto i = 0; i < n; i++)  // Iterate through elements
        v[i] = 10*fabsf64x(rand()%100+11); // Set element to random value in [10, 100]
    return v;  // Return the generated vector
}

int main(int argc, char** argv){
	const int N = 5;
    
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

    // Uncomment this code to test the example; The answer is ~=(1,1.92,-1,3,-0.8)
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

	simfor::vec res_vec1 = simfor::SteepestDescentSolverOmp(myMat, myVec, 20*N, 1e-6);
	simfor::vec res_vec2 = simfor::ConjugateGradientSolverOmp(myMat, myVec, 20*N, 1e-6);

    std::cout << "Answer(steepestDescentSolver): " << [&res_vec1](){ for (auto &&i : res_vec1){std::cout << i << " ";}; return "\n";}();

    std::cout << "Answer(ConjugateGradientSolver): " << [&res_vec2](){ for (auto &&i : res_vec2){std::cout << i << " ";}; return "\n";}();

	return 0;
}
#include "simfor/LUdecompOmp.hpp"

/**
 * @brief Generate a matrix with random values
 * 
 * @param n size of the matrix
 * @return simfor::matr n x n matrix with random values
 */
simfor::matr genMatNNB(int n)
{
    simfor::matr m(n, n); // create a n x n matrix
    for (auto i = 0; i < n; i++) // iterate over rows
    {
        for (auto j = 0; j < n; j++) // iterate over columns
        {
            if (i == j) // if it is a diagonal element set a random value between 10 and 20
            {
                m(i, j) = 10 * fabsf64x(rand() % 10 + 11);
            }
            else // otherwise set a random value between 0 and 9
            {
                m(i, j) = rand() % 10;
            }
        }
    }
    return m;
}

/**
 * @brief Generate a vector with random values
 * 
 * @param n Size of the vector
 * @return simfor::vec Vector of length n with random values between 10 and 20
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); // Create a vector with length n
    for(auto i = 0; i < n; /* increment i */) {
        v[i] = 10*fabsf64x(rand()%10+11); // Set the current element to a random value between 10 and 20
        i++; // Increment i
    }
    return v;
}

int main(int argc, char** argv){
	const auto N = 5;

    simfor::vec P(N+1), X(N);
    
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
    
    auto status = simfor::LUPDecomposeOmp(myMat, size_t(N), 1e-6, P);
    simfor::LUPSolveOmp(myMat, P, myVec, size_t(N), X);

    std::cout << "Answer: " << [&X](){ for (auto &&i : X){std::cout << i << " ";}; return "\n";}();

	return 0;
}
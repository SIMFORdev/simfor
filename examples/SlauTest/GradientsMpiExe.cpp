#include "simfor/GradientsMpi.hpp"

/**
 * @brief Generate a random matrix with n rows and n columns
 * @param n The number of rows and columns
 * @return A matrix with a diagonal of random large elements and the rest of random small elements
 */
simfor::matr genMatNNB(int n){
        simfor::matr m(n, n); /* create a matrix with n rows and n columns */
        /* fill the matrix with random values */
        for(auto i=0;i<n;i++){ /* iterate over rows */
            for(auto j=0;j<n;j++){ /* iterate over columns */
                if (i==j) /* if the current element is on the diagonal */
                {
                    m(i,j) = 10*fabsf64x(rand()%100+11); /* set a random large value */
                }else{
                    m(i,j) = rand()%10; /* set a random small value */
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector of length n with random values
 * @param n The length of the vector
 * @return A vector with random values between 10 and 110
 */
simfor::vec genVecN(int n){
    simfor::vec v(n); /* create a vector of length n */

    /* fill the vector with random values */
    for(auto i = 0; i < n; /* increment i and do the following */)
        v[i++] = 10*fabsf64x(rand()%10+11); /* set i'th element to a random value between 10 and 110 */

    return v;
}

int main(int argc, char** argv){
    mpi::environment env (argc, argv);
    mpi::communicator world;

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

	simfor::vec res_vec1 = simfor::SteepestDescentSolverMpi(myMat, myVec, 20*N, 1e-6);
	simfor::vec res_vec2 = simfor::ConjugateGradientSolverMpi(myMat, myVec, 20*N, 1e-6);

    std::cout << "Answer(steepestDescentSolver): " << [&res_vec1](){ for (auto &&i : res_vec1){std::cout << i << " ";}; return "\n";}();

    std::cout << "Answer(ConjugateGradientSolver): " << [&res_vec2](){ for (auto &&i : res_vec2){std::cout << i << " ";}; return "\n";}();

	return 0;
}
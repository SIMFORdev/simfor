#include "simfor/SimpleIterMpi.hpp"

/**
 * @brief Generate a random square matrix with random non-zero diagonal elements
 * @param n Size of the matrix
 * @return A square matrix of size n, with diagonal elements randomized
 */
simfor::matr genMatNNB(int n){
        /* Generate a random square matrix with random non-zero diagonal elements */
        simfor::matr m(n, n); /* Create a matrix of size n x n */
        for(auto i=0;i<n;i++){ /* Iterate over rows */
            for(auto j=0;j<n;j++){ /* Iterate over columns */
                if (i==j) /* If we're on the diagonal */
                {
                    m(i,j) = fabsf64x(rand()%10+11); /* Set diagonal element to a random value */
                }else{
                    m(i,j) = rand()%10; /* Otherwise set element to a random value */
                }
            }
        }
        return m; /* Return the generated matrix */
}

/**
 * @brief Generate a vector with n elements
 * @param n Number of elements in the vector
 * @return A vector with n elements, with random values in [10, 100]
 */
simfor::vec genVecN(int n)
{
    simfor::vec v(n); // Create a vector of size n
    for (auto i = 0; i < n; ) // Iterate over the vector
    {
        v[i++] = 10*fabsf64x(rand()%10+11); // Set each element to a random value between [10, 100]
    }
    return v; // Return the generated vector
}

int main(int argc, char** argv){
	const auto N = 5;
    mpi::environment env (argc, argv);
    mpi::communicator world;
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

    simfor::vec resVec = simfor::SimpleIterMpi(myMat, myVec, N);

    std::cout << "Answer: " << [&resVec](){ for (auto &&i : resVec){std::cout << i << " ";}; return "\n";}();

	return 0;
}
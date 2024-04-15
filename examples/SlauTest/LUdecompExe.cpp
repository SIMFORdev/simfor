#include "simfor/LUdecomp.hpp"

/**
 * @brief Generate a matrix with a diagonal of random numbers between 10 and 20
 * and the rest of the numbers random between 0 and 9.
 *
 * @param n The size of the matrix
 *
 * @return A matrix with the described properties
 */
simfor::matr genMatNNB(int n){
        simfor::matr m(n, n); // Create an n x n matrix
        for(int i=0;i<n;i++){ // Iterate over rows
            for(int j=0;j<n;j++){ // Iterate over columns
                if (i==j)
                {
                    m(i,j) = 10*fabsf64x(rand()%10+11); // Set diagonal elements to random numbers between 10 and 20
                }else{
                    m(i,j) = rand()%10; // Set other elements to random numbers between 0 and 9
                }
            }
        }
        return m;
}

/**
 * @brief Generate a vector with n elements
 *
 * The vector is filled with random numbers between 10 and 20
 *
 * @param n Number of elements in the vector
 *
 * @return A vector with n elements, each element is a random number between 10 and 20
 */
simfor::vec genVecN(int n){
    simfor::vec v(n);  /* create a vector with n elements */
    
    /* fill the vector with random numbers between 10 and 20 */
    for(auto i = 0; i < n; /* increment i after the loop body */)
    {
        v[i] = 10*fabsf64x(rand()%10+11);  /* set the current element to a random number between 10 and 20 */
        ++i;  /* increment i before the loop body */
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
    
    auto status = simfor::LUPDecompose(myMat, size_t(N), 1e-6, P);
    simfor::LUPSolve(myMat, P, myVec, size_t(N), X);

    std::cout << "Answer: " << [&X](){ for (auto &&i : X){std::cout << i << " ";}; return "\n";}();

	return 0;
}
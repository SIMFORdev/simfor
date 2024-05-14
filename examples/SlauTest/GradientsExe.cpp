#include "simfor/Gradients.hpp"

/**
 * @brief Generate matrix with n x n size filled with numbers
 * 
 * @param n - matrix size
 * @return simfor::matr - generated matrix with numbers
 */
simfor::matr genMatNNB(int n){
        simfor::matr m(n, n); // initialize matrix with n x n size
        for(auto i=0;i<n;i++){ // iterate through rows
            for(auto j=0;j<n;j++){ // iterate through columns
                if (i==j) // if element is on diagonal
                {
                    m(i,j) = 10*fabsf64x(rand()%10+11); // set element to 10 * abs(random number) + 11
                }else{
                    m(i,j) = rand()%10; // set element to random number
                }
            }
        }
        return m;
} // end of function

/**
 * @brief Generate vector with n elements filled with numbers
 * 
 * @param n - vector size
 * @return simfor::vec - generated vector with numbers
 */
simfor::vec genVecN(int n){
    /* Initialize vector with n elements */
    simfor::vec v(n);
    /* Iterate through the vector */
    for(auto i = 0; i < n; /* increment inside loop */)
    {
        /* Set i-th element to 10 * abs(random number) + 11 */
        v[i++] = 10*fabsf64x(rand()%10+11);
    }
    /* Return generated vector */
    return v;
} // end of function genVecN

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

	simfor::vec res_vec1 = simfor::SteepestDescentSolver(myMat, myVec, 20*N, 1e-6);
	simfor::vec res_vec2 = simfor::ConjugateGradientSolver(myMat, myVec, 20*N, 1e-6);

    std::cout << "Answer(steepestDescentSolver): " << [&res_vec1](){ for (auto &&i : res_vec1){std::cout << i << " ";}; return "\n";}();

    std::cout << "Answer(ConjugateGradientSolver): " << [&res_vec2](){ for (auto &&i : res_vec2){std::cout << i << " ";}; return "\n";}();

	return 0;
}
#include "simfor/Gradients.hpp"

// simfor::matr genMatNNB(int n){
//         simfor::matr m(n, n);
//         for(auto i=0;i<n;i++){
//             for(auto j=0;j<n+1;j++){
//                 if (i==j)
//                 {
//                     m(i,j) = fabsf64x(rand()%10+11);
//                 }else{
//                     m(i,j) = rand()%10;
//                 }
//             }
//         }
//         return m;
// }

// simfor::vec genVecN(int n){
//     simfor::vec v(n);
//     for(auto i = 0; i < n; v[i++] = 10*rand()%10); 
//     return v;
// }

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
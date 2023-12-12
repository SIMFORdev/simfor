#include "simfor/SimpleIterOmp.hpp"

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
    simfor::matr myMat(5, 5);
    simfor::vec myVec(5);

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
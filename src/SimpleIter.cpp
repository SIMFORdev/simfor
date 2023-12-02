#include "simfor/SimpleIter.hpp"

namespace simfor{

    vec SimpleIter(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N);
        long double eps = 1e-9;
        size_t iter{};
        long double sum{};
        
        for (auto i = 0; i < N; i++){
            x[i] = vecB[i] / mat(i,i);
        }

        do{
            // flag = 1;
            for(auto i = 0; i < N; i++){
                x_n[i] = vecB[i] / mat(i,i);
                for(auto j = 0; j < N; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n[i] -= mat(i,j) / mat(i,i) * x[j];
                    };
                }
            }

            int flag = 1;
            for(auto i = 0; i < N-1; i++){
                if (std::abs(x_n[i] - x[i]) > eps){
                    flag = 0;
                    break;
                }
            }

            for (auto i = 0; i < N; i++){
                x[i] = x_n[i];
            }

            if (flag) break;

        } while (1);

        return x;
    }


}
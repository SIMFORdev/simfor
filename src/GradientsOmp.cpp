#include <cmath>           //
#include <numeric>         //inner_product - скалярное произведение векторов       

#include "omp.h"
#include "simfor/GradientsOmp.hpp"

namespace simfor{

    //метод сопряжённых градиентов
    vec ConjugateGradientSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps){
        auto pcnt = omp_get_max_threads();
        omp_set_num_threads(pcnt);
        auto n = A.size1(); //возвращает кол-во строк матрицы
        vec X(n, 0.0); //вектор решения
        vec R = B;
        vec P = R;
        size_t k{};  //кол-во текущих итераций

        while ((k < iter_num) && (vectorNorm(R) > eps)){
            if (k == 0)
                P = R;
            else{
                auto gamma = -(innerProduct(P, matrixMulVector(A, R)))/(innerProduct(P, matrixMulVector(A, P)));
                P = vectorCombination(1.0, R, gamma, P);
            }
            auto alpha = innerProduct(P, R) / innerProduct(P, matrixMulVector(A, P));
            X = vectorCombination(1.0, X, alpha, P);                                        // След. ожидаемое значение
            R = vectorCombination(1.0, R, -alpha, matrixMulVector(A, P)); 
            k++;
            if (k >= iter_num){
                std::cout << "(сопр. градиент)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
                break;
            }
        }

        return X;
    }

    vec SteepestDescentSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps){
        auto pcnt = omp_get_max_threads();
        omp_set_num_threads(pcnt);
        // std::cout << "Threads in work" << pcnt << '\n';
        auto n = A.size1(); size_t k{};
        vec X(n, 0.0);
        vec R = B;

        while ((k < iter_num) && (vectorNorm(R) > eps)){
            vec P = R;
            vec Q = matrixMulVector(A, P);
            double alpha, a1, a2;

            a1 = innerProduct(P, R);

            a2 = innerProduct(P, Q);

            alpha = a1/a2;
            X = vectorCombination(1.0, X, alpha, P);            // След. ожидаемое значение
                        
            R = vectorCombination(1.0, R, -alpha, Q);          // Остаток

            k++;
            if ((k >= iter_num)){
               std::cout << "(наиск. спуск)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
                break;
            }
        }
        return X;
    }

    vec matrixMulVector(const matr &A, const vec &V)     // умножение матрицы на вектор
    {
        auto n = A.size1();
        vec C(n), tmp(n);
        size_t i{}, j{};
        #pragma omp taskloop default(shared) private(i)
        for (i = 0; i < n; i++){
            #pragma omp parallel for private(j)
            for (j = 0; j < n; j++){
                tmp[j] = A(i, j);
            }
            C[i] = innerProduct(tmp, V);
        }
        return C;
    }

    double vectorNorm(const vec &V){
        return sqrt(innerProduct(V, V));
    }

    double innerProduct(const vec &U, const vec &V){
        // return inner_product(U.begin(), U.end(), V.begin(), 0.0);
        double res{};
        size_t i{};

        #pragma omp parallel for private(i) reduction(+:res)
        for (i = 0; i < U.size(); i++){
            res += U[i] * V[i];
        }
        
        return res;
    }

    vec vectorCombination(double a, const vec &U, double b, const vec &V)        // Linear combination of vectors
    {
        auto n = U.size();
        vec W(n);
        size_t i{};

        #pragma omp parallel for private(i) 
        for (i = 0; i < n; i++ ) W[i] = a * U[i] + b * V[i];

        return W;
    }

}
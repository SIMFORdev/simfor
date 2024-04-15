#include <cmath>           //
#include <numeric>         //inner_product - скалярное произведение векторов       

#include "omp.h"
#include "simfor/GradientsOmp.hpp"

namespace simfor{

    /**
     * @brief Метод сопряженных градиентов
     *
     * Решение системы линейных уравнений Ax=b посредством метода сопряженных градиентов.
     *
     * @param A коэффициенты матрицы А
     * @param B вектор правой части
     * @param iter_num максимальное количество итераций
     * @param eps точность решения
     * @return вектор решения
     */
    vec ConjugateGradientSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps){
        auto pcnt = omp_get_max_threads();
        omp_set_num_threads(pcnt);
        auto n = A.size1(); //возвращает кол-во строк матрицы
        vec X(n, 0.0); //вектор решения
        vec R = B;
        vec P = R;
        size_t k{};  //кол-во текущих итераций

        while ((k < iter_num) && (vectorNorm(R) > eps)){
            if (k == 0)   //если первая итерация
                P = R;    //P = R
            else{
                //gamma = -(R^T A R) / (P^T A P)
                auto gamma = -(innerProduct(P, matrixMulVector(A, R)))/(innerProduct(P, matrixMulVector(A, P)));
                P = vectorCombination(1.0, R, gamma, P);
            }
            //alpha = R^T A P / P^T A P
            auto alpha = innerProduct(P, R) / innerProduct(P, matrixMulVector(A, P));
            X = vectorCombination(1.0, X, alpha, P); // X = X + alpha P
            R = vectorCombination(1.0, R, -alpha, matrixMulVector(A, P)); // R = R - alpha A P
            k++;
            if (k >= iter_num){
                std::cout << "(сопр. градиент)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
                break;
            }
        }

        return X;
    }


    /**
     * Метод наискорейшего спуска
     *
     * Решение системы линейных уравнений Ax=b посредством метода наискорейшего спуска.
     *
     * @param A коэффициенты матрицы А
     * @param B вектор правой части
     * @param iter_num максимальное количество итераций
     * @param eps требуемая точность
     * @return вектор решения
     */
    vec SteepestDescentSolverOmp(const matr &A, const vec &B, size_t iter_num, double eps){
        auto pcnt = omp_get_max_threads();
        omp_set_num_threads(pcnt);
        auto n = A.size1(); size_t k{};
        vec X(n, 0.0); // вектор решения
        vec R = B;     // вектор ошибки

        while ((k < iter_num) && (vectorNorm(R) > eps)){
            vec P = R; // вектор итерации
            vec Q = matrixMulVector(A, P); //A P
            double alpha, a1, a2;         // коэффициенты<|context_not_required|>

            a1 = innerProduct(P, R);  // <P, R>

            a2 = innerProduct(P, Q);  // <P, A P>

            alpha = a1/a2; // коэффициент итерации
            X = vectorCombination(1.0, X, alpha, P); // X = X + alpha P
                                                      // След. ожидаемое значение
            R = vectorCombination(1.0, R, -alpha, Q); // R = R - alpha A P
                                                      // Остаток

            k++;
            if ((k >= iter_num)){
                std::cout << "(наиск. спуск)Решение не может быть достигнуто при точности: " << eps << " за "
                          << iter_num << " итераций\n";
                break;
            }
        }
        return X;
    }


    vec matrixMulVector(const matr &A, const vec &V)     // умножение матрицы на вектор
    {
        auto n = A.size1();
        vec C(n);
        const double* v = &V[0];

        #pragma omp parallel for default(none) shared(A, v, C) firstprivate(n)
        for (size_t i = 0; i < n; ++i) {
            double x = 0;
            const double* a_i = &A(i, 0);
            for (size_t j = 0; j < n; ++j) {
                x += a_i[j] * v[j];
            }
            C[i] = x;
        }
        return C;
    }


    double vectorNorm(const vec &V){
        return sqrt(innerProduct(V, V));
    }

    double innerProduct(const vec &U, const vec &V){
        const auto n = U.size();
        double res{};

        const double *pu = &U[0];
        const double *pv = &V[0];

        #pragma omp simd reduction(+:res) aligned(pu,pv:sizeof(double))
        for (std::size_t i = 0; i < n; ++i) res += pu[i] * pv[i];

        return res;
    }


    vec vectorCombination(double a, const vec &U, double b, const vec &V)        // Linear combination of vectors
    {
        auto n = U.size();
        vec W(n);

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++ ) W[i] = a * U[i] + b * V[i];

        return W;
    }


}
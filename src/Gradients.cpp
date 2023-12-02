#include <cmath>           //
#include <numeric>         //inner_product - скалярное произведение векторов
#include "simfor/Gradients.hpp"

namespace simfor{

   vec ConjugateGradientSolver(const matr &A, const vec &B, size_t iter_num, double eps){
      //метод сопряжённых градиентов
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
         double alpha = innerProduct(P, R) / innerProduct(P, matrixMulVector(A, P));
         X = vectorCombination(1.0, X, alpha, P);            // След. ожидаемое значение
         R = vectorCombination(1.0, R, -alpha, matrixMulVector(A, P)); 
         k++;
         if ((k >= iter_num)){
            std::cout << "(сопр. градиент)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
            break;
         }
      }
      return X;
   }

   vec SteepestDescentSolver(const matr &A, const vec &B, size_t iter_num, double eps){
      auto n = A.size1();
      vec X(n, 0.0);

      vec R = B;
      size_t k{};

      while ((k < iter_num) && (vectorNorm(R) > eps)){
         vec P = R;
         vec Q = matrixMulVector(A, P);
         auto alpha = innerProduct(P, R) / innerProduct(P, Q);
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
      vec C(n);
      vec tmp(n);
      for (auto i = 0; i < n; i++){
         for (auto j = 0; j < n; j++){
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
      return inner_product(U.begin(), U.end(), V.begin(), 0.0);
   }

   vec vectorCombination(double a, const vec &U, double b, const vec &V)        // Linear combination of vectors
   {
      auto n = U.size();
      vec W(n);
      for (auto j = 0; j < n; j++ ) W[j] = a * U[j] + b * V[j];
      return W;
   }

}
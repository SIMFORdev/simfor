#include <cmath>           //
#include <numeric>         //inner_product - скалярное произведение векторов
#include <simfor/Gradients.hpp>

namespace simfor{

   /**
    * Метод сопряженных градиентов
    * @param A - матрица коэффициентов
    * @param B - правая часть системы уравнений
    * @param iter_num - максимальное количество итераций
    * @param eps - точность, до которой нужно достичь
    * @return вектор решения
    */
   vec ConjugateGradientSolver(const matr &A, const vec &B, size_t iter_num, double eps){
      auto n = A.size1(); //возвращает кол-во строк матрицы
      vec X(n, 0.0); //вектор решения

      vec R = B; // вектор ошибки
      vec P = R; // вектор направления поиска
      size_t k{};  //кол-во текущих итераций

      while ((k < iter_num) && (vectorNorm(R) > eps)){
         if (k == 0)
            P = R; // начальное направление поиска
         else{
            /**
             * коэффициент сопряженности, равный отношению косинусов
             * угла между вектором направления поиска на предыдущем
             * шаге и вектором, полученным из умножения матрицы на вектор
             * ошибки на предыдущем шаге
             */
            auto gamma = -(innerProduct(P, matrixMulVector(A, R)))/(innerProduct(P, matrixMulVector(A, P)));
            P = vectorCombination(1.0, R, gamma, P);
         }
         /**
          * коэффициент на шаге итерации, равный отношению косинуса
          * угла между вектором направления поиска на текущем
          * шаге и вектором, полученным из умножения матрицы на вектор
          * направления поиска на текущем шаге
          */
         double alpha = innerProduct(P, R) / innerProduct(P, matrixMulVector(A, P));
         X = vectorCombination(1.0, X, alpha, P);            // След. ожидаемое значение
         R = vectorCombination(1.0, R, -alpha, matrixMulVector(A, P)); // Остаток
         k++;
         if ((k >= iter_num)){
            std::cout << "(сопр. градиент)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
            break;
         }
      }
      return X;
   }


   /**
    * Метод наискорейшего спуска
    *
    * @param A матрица коэффициентов
    * @param B вектор правой части
    * @param iter_num максимальное количество итераций
    * @param eps точность
    * @return вектор решения
    */
   vec SteepestDescentSolver(const matr &A, const vec &B, size_t iter_num, double eps){
      auto n = A.size1(); // размерность матрицы
      vec X(n, 0.0); // вектор решения

      vec R = B; // вектор ошибки
      size_t k{};  // кол-во текущих итераций

      while ((k < iter_num) && (vectorNorm(R) > eps)){
         vec P = R; // вектор направления поиска
         vec Q = matrixMulVector(A, P); // умножение матрицы на вектор направления поиска
         auto alpha = innerProduct(P, R) / innerProduct(P, Q); // коэффициент на шаге итерации
         X = vectorCombination(1.0, X, alpha, P); // След. ожидаемое значение
         R = vectorCombination(1.0, R, -alpha, Q); // Остаток
         k++;
         if ((k >= iter_num)){
            std::cout << "(наиск. спуск)Решение не может быть достигнуто при точности: "
                      << eps << " за " << iter_num << " итераций\n";
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
            tmp(j) = A(i, j);
         }
         C(i) = innerProduct(tmp, V);
      }
      return C;
   }

   double vectorNorm(const vec &V){
      return sqrt(innerProduct(V, V));
   }

   double innerProduct(const vec &U, const vec &V){
      return inner_product(U.begin(), U.end(), V.begin(), 0.0);
   }

   /**
    * Линейная комбинация векторов
    *
    * @param a коэффициент вектора U
    * @param U вектор U
    * @param b коэффициент вектора V
    * @param V вектор V
    * @return вектор, являющийся линейной комбинацией U и V
    */
   vec vectorCombination(double a, const vec &U, double b, const vec &V){
      auto n = U.size();
      vec W(n);

      // формула линейной комбинации векторов:
      // W = a * U + b * V
      std::transform(U.begin(), U.end(), V.begin(), W.begin(),
                     [a, b](double ui, double vi){ return a * ui + b * vi; });

      return W;
   }


}
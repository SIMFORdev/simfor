#include <cmath>           //
#include <numeric>         //inner_product - скалярное произведение векторов
#include "simfor/GradientsMpi.hpp"
#include <iostream>

namespace simfor{

   vec ConjugateGradientSolverMpi(const matr &A, const vec &B, size_t iter_num, double eps){
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
            auto gamma = -(innerProduct(P, matrixMulVectorMpi(A, R)))/(innerProduct(P, matrixMulVectorMpi(A, P)));
            P = vectorCombinationMpi(1.0, R, gamma, P);
         }
         double alpha = innerProduct(P, R) / innerProduct(P, matrixMulVectorMpi(A, P));
         X = vectorCombinationMpi(1.0, X, alpha, P);            // След. ожидаемое значение
         R = vectorCombinationMpi(1.0, R, -alpha, matrixMulVectorMpi(A, P)); 
         k++;
         if ((k >= iter_num)){
            std::cout << "(сопр. градиент)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
            break;
         }
      }
      return X;
   }

   vec SteepestDescentSolverMpi(const matr &A, const vec &B, size_t iter_num, double eps){
        auto n = A.size1();
        vec X(n, 0.0);

        vec R = B;
        size_t k{};

        mpi::communicator world;
        int p = world.size();
        int r = world.rank();

        int cnt = floor(n / p);
        int from = r * cnt;
        int to = n;
        if ( r != p-1 ) to = from + cnt;

        while ((k < iter_num) && (vectorNorm(R) > eps)){
            world.barrier();
            double alpha, a1, a2;
            vec P = R;
            vec Q = matrixMulVectorMpi(A, P);

            a1 = innerProduct(P, R);
            a2 = innerProduct(P, Q);
            alpha = a1/a2;
            X = vectorCombinationMpi(1.0, X, alpha, P);            // След. ожидаемое значение
            R = vectorCombinationMpi(1.0, R, -alpha, Q);          // Остаток
            k++;
            if ((k >= iter_num)){
                std::cout << "(наиск. спуск)Решение не может быть достигнуто при точности: " << eps << " за " << iter_num << " итераций\n";
                break;
            }
            world.barrier();
        }
        return X;
   }

   vec matrixMulVectorMpi(const matr &A, const vec &V)     // умножение матрицы на вектор
   {
        auto n = A.size1();
        vec C(n, 0.0);
        vec tmp(n, 0.0);

        mpi::communicator world;
        int p = world.size();
        int r = world.rank();

        int cnt = floor(n / p);
        int from = r * cnt;
        int to = n;
        if ( r != p-1 ) to = from + cnt;

        for (auto i = from; i < to; i++){
            for (auto j = 0; j < n; j++){
                tmp(j) = A(i, j);
            }
            C(i) = innerProduct(tmp, V);
        }

        world.barrier();
        if (r != 0){
            world.send( 0, 1, from );
            world.send( 0, 2, to );
            world.send(0, 3, &C(from), ( to - from));
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &C(from), ( to - from));
            }
        }
        broadcast(world, &C(0), n,  0);
        return C;
   }

   double vectorNorm(const vec &V){
      return sqrt(innerProduct(V, V));
   }

   double innerProduct(const vec &U, const vec &V){
        return inner_product(U.begin(), U.end(), V.begin(), 0.0);
        // double res{0.0}, tmp{0.0};
        // auto n = U.size();

        // mpi::communicator world;
        // int p = world.size();
        // int r = world.rank();
        // int cnt = floor(n / p);
        // int from = r * cnt;
        // int to = n;
        // if ( r != p-1 ) to = from + cnt;

        // for (int i = from; i < to; ++i){
        //     tmp += U(i) * V(i);
        // }

            // reduce(world, tmp, res, std::plus<float>(), 0);
        // broadcast(world, &res, 1, 0);
        // world.barrier();
        // return res;
   }

   vec vectorCombinationMpi(double a, const vec &U, double b, const vec &V)        // Linear combination of vectors
   {
        auto n = U.size();
        vec W(n, 0.0);

        mpi::communicator world;
        int p = world.size();
        int r = world.rank();
        int cnt = floor(n / p);
        int from = r * cnt;
        int to = n;
        if ( r != p-1 ) to = from + cnt;

        for (auto j = from; j < to; j++ ) W(j) = a * U(j) + b * V(j);

        world.barrier();
        if (r != 0){
            world.send( 0, 1, from );
            world.send( 0, 2, to );
            world.send(0, 3, &W(from), ( to - from));
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &W(from), ( to - from));
            }
        }
        broadcast(world, &W(0), n,  0);

        return W;
   }

}
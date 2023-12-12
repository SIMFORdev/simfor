#include <iostream>
#include "simfor/GaussOmp.hpp"

namespace simfor{

	vec GaussianEliminationOmp(matr &mat, int N){
		/*Сведение матрицы к треугольной*/
		int singular_flag = ForwardElimOmp(mat, N);
		vec res;
		// printf("After forward elim:\n");
		// print(mat, N);
		/* Если матрица вырожд */
		if (singular_flag != -1){
			printf("Singular Matrix.\n");
			/* если правая часть уравнения, 
			соответствующая нулевой строке, 
			равна 0, * система имеет бесконечно много решений, 
			в противном случае система несовместна*/
			if (mat(singular_flag, N))
				printf("Inconsistent System.");
			else
				printf("May have infinitely many solutions.");
			return res;
		}
		/* get solution to system and print it using
		backward substitution */
		res = BackSubOmp(mat, N);
		return res;
	}

	void SwapRow(matr &mat, int i, int j, int N){
		for (auto k = 0; k <= N; k++){
			double temp = mat(i, k);
			mat(i, k) = mat(j, k);
			mat(i, k) = temp;
		}
	}

	int ForwardElimOmp(matr &mat, int N){
		int i{}, j{}, k{}, i_max{}, v_max{}, flag{};
		double f{};

		// #pragma omp parallel for default(shared) private(k, i_max, v_max, i, j, f)
		#pragma omp taskloop default(shared) private(k, i_max, v_max, i, j, f) // num_tasks(N)
		for (k=0; k<N; k++){
			// Инициализируем максимальное значение и индекс
			i_max = k;	//index
			v_max = mat(i_max, k);	//value

			/* находим наибольш знач для свапа, если таковые имеются */
			// #pragma omp parallel for default(shared) private(i)
			// #pragma omp single nowait
			#pragma omp simd order(concurrent)// nontemporal(mat)
			for (i = k+1; i < N; i++)
				if (std::abs(mat(i, k)) > v_max)
					v_max = mat(i, k), i_max = i;

			/* если элемент главной диагонали близок нулю*, это означает,
			что матрица вырожд., что приведет к потере точности. */
			if (std::abs(mat(k, i_max)) < 1e-9)
				// return k;
				exit(-1);

			/* Поменять местами строку с наибольшим значением на текущую строку */
			if (i_max != k)
				SwapRow(mat, k, i_max, N);

			// #pragma omp taskloop default(shared) private(i,j,f)
			#pragma omp parallel for default(shared) private(i,j,f)
			for (i=k+1; i<N; i++){
				/* коэффициент f*/
				f = mat(i, k)/mat(k, k);

				/* вычесть f-к кратное соответствующего k-го элемента строки*/
				// #pragma omp prallel for default(shared) private(j)
				#pragma omp simd order(concurrent)
				for(j = k+1; j<=N; j++)
					mat(i, j) -= mat(k, j)*f;

				/* нижний трегуольник заполняем нулями*/
				mat(i, k) = 0;
			}
		}
		return -1;
	}

	vec BackSubOmp(matr &mat, int N){
		vec x(N);
		int i{},j{};
		long double tmp{};

		/* Начнём расчет с последнего уравнения до первого */
		#pragma omp taskloop default(shared) private(i,j,tmp)
		for (i = N-1; i >= 0; i--){
			/* начнём с правой стороны ур-ия */
			tmp = mat(i, N);
			/* Инициализируйте j значением i+1, так как матрица является верхнетреугольной.*/
			#pragma omp simd order(concurrent) private(j) reduction(-:tmp) 
			for (j=i+1; j<N; j++){
				/* вычтим всем значения
				* кроме коэффициента при переменной
				* значение которого рассчитывается */
				tmp -= mat(i, j)*x(j);
			}
			/* разделим правую часть ур-ия на коэффициент при вычисляемом неизвестном */
			x(i) = tmp/mat(i, i);
		}
		return x;
	}
}
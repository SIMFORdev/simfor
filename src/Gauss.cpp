#include <iostream>
#include "simfor/Gauss.hpp"

namespace simfor{

	vec GaussianElimination(matr &mat, int N){
		/*Сведение матрицы к треугольной*/
		int singular_flag = ForwardElim(mat, N);
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
		res = BackSub(mat, N);
		return res;
	}

	void SwapRow(matr &mat, int i, int j, int N){
		for (auto k = 0; k <= N; k++){
			double temp = mat(i, k);
			mat(i, k) = mat(j, k);
			mat(i, k) = temp;
		}
	}



	int ForwardElim(matr &mat, int N){
		for (auto k=0; k<N; k++){
			// Инициализируем максимальное значение и индекс
			auto i_max = k;	//index
			auto v_max = mat(i_max, k);	//value

			/* находим наибольш знач для свапа, если таковые имеются */
			for (auto i = k+1; i < N; i++)
				if (std::abs(mat(i, k)) > v_max)
					v_max = mat(i, k), i_max = i;

			/* если элемент главной диагонали близок нулю*, это означает,
			что матрица вырожд., что приведет к потере точности. */
			if (std::abs(mat(k, i_max)) < 1e-9)
				return k;

			/* Поменять местами строку с наибольшим значением на текущую строку */
			if (i_max != k)
				SwapRow(mat, k, i_max, N);

			for (auto i=k+1; i<N; i++){
				/* коэффициент f*/
				double f = mat(i, k)/mat(k, k);

				/* вычесть f-к кратное соответствующего k-го элемента строки*/
				for(auto j=k+1; j<=N; j++)
					mat(i, j) -= mat(k, j)*f;

				/* нижний трегуольник заполняем нулями*/
				mat(i, k) = 0;
			}
		}
		return -1;
	}

	vec BackSub(matr &mat, int N){
		vec x(N);
		/* Начнёс расчет с последнего уравнения до первого */
		for (auto i = N-1; i >= 0; i--){
			/* начнём с правой стороны ур-ия */
			x[i] = mat(i, N);
			/* Инициализируйте j значением i+1, так как матрица является верхнетреугольной.*/
			for (auto j=i+1; j<N; j++){
				/* вычтим всем значения
				* кроме коэффициента при переменной
				* значение которого рассчитывается */
				x[i] -= mat(i, j)*x[j];
			}
			/* разделим правую часть ур-ия на коэффициент при вычисляемом неизвестном */
			x[i] /= mat(i, i);
		}
		return x;
	}

}
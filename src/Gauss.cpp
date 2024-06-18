#include <iostream>
#include <algorithm>
#include <simfor/Gauss.hpp>

namespace simfor{

	/**
	 * Performs Gaussian Elimination on the given matrix and right-hand side vector
	 * @param[in,out] mat The matrix to eliminate, modified in place.
	 * @param[in] N The size of the square matrix.
	 * @return The solution vector to the system of linear equations.
	 * 
	 * Gaussian Elimination reduces the given matrix (A) and right-hand side vector (b)
	 * to row-echelon form (REQ). The solution to the system of linear equations is
	 * then obtained by performing backward substitution on the modified matrix (A) and
	 * vector (b).
	 * 
	 * If the matrix is singular, the function prints out a message indicating
	 * whether the system is inconsistent or may have infinitely many solutions. In
	 * either case, an empty vector is returned.
	 */
	vec GaussianElimination(matr &mat, int N){
		// Perform Gaussian Elimination
		int singular_flag = ForwardElim(mat, N);
		// If matrix is singular, handle accordingly
		if (singular_flag != -1){
			printf("Singular Matrix.\n");
			// If right-hand side is zero, infinitely many solutions
			if (mat(singular_flag, N))
				printf("Inconsistent System.\n");
			// Otherwise may have infinitely many solutions
			else
				printf("May have infinitely many solutions.\n");
			return vec();
		}
		
		// Solve the system of linear equations using backward substitution
		vec res = BackSub(mat, N);
		return res;
	}


	void SwapRow(matr &mat, int i, int j, int N){
		std::swap_ranges(mat.begin1() + i * (N + 1), mat.begin1() + (i + 1) * (N + 1), mat.begin1() + j * (N + 1));
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
			x(i) = mat(i, N);
			/* Инициализируйте j значением i+1, так как матрица является верхнетреугольной.*/
			for (auto j=i+1; j<N; j++){
				/* вычтим всем значения
				* кроме коэффициента при переменной
				* значение которого рассчитывается */
				x(i) -= mat(i, j)*x(j);
			}
			/* разделим правую часть ур-ия на коэффициент при вычисляемом неизвестном */
			x(i) /= mat(i, i);
		}
		return x;
	}

}
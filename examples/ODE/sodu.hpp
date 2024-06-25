
#ifndef SIMFOR_SODU_HPP
#define SIMFOR_SODU_HPP

#include <omp.h>
#include <iostream>
#include "multMatrixVector.hpp"


namespace simfor
{
//коэффаценты соотв. методов
double adams_bashford_koeff( int );
double adams_moulton_koeff( int );

//метотоды Адамса-Милна
matrix<double> AdMiln_system_solve_matrix ( double, int, vector<double>, matrix<double> );
matrix<double> AdMiln_system_solve_matrix_omp ( double, int, vector<double>, matrix<double> );
matrix<double> AdMiln_system_solve_matrix_mpi ( double, int, vector<double>, matrix<double> );

//метотоды Адамса-БашфортаW
matrix<double> adams5_system_solve_matrix ( double, int, vector<double>, matrix<double> );
matrix<double> adams5_system_solve_matrix_omp ( double, int, vector<double>, matrix<double> );
matrix<double> adams5_system_solve_matrix_mpi ( double, int, vector<double>, matrix<double> );

//метотоды Адамса-Моултона

matrix<double> AdMltn_system_solve_matrix ( double, int, vector<double>, matrix<double> );
matrix<double> AdMltn_system_solve_matrix_omp ( double, int, vector<double>, matrix<double> );
matrix<double> AdMltn_system_solve_matrix_mpi ( double, int, vector<double>, matrix<double> );

//eiler
matrix<double> eiler_system_solve_matrix ( double, int, vector<double>, matrix<double>);
matrix<double> eiler_system_solve_matrix_omp ( double, int, vector<double>, matrix<double>);
matrix<double> eiler_system_solve_matrix_mpi ( double, int, vector<double>, matrix<double>);

//rk(Runge-Kutta)
matrix<double> rk_system_solve_matrix ( double, int, vector<double>, matrix<double>);
matrix<double> rk_system_solve_matrix_omp ( double, int, vector<double>, matrix<double>);
matrix<double> rk_system_solve_matrix_mpi ( double, int, vector<double>, matrix<double>);
}

#endif

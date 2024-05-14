#ifndef SIMFOR_ODU_HPP
#define SIMFOR_ODU_HPP

#include <omp.h>
#include <boost/mpi.hpp>
#include "internal/types.hpp"
#include "elementary.hpp"
#include "simfor/multMatrVec.hpp"
#include <iostream>



namespace simfor
{
namespace mpi = boost::mpi;
//Coeffisients for Adams methods
float adams_bashford_koeff( int );
float adams_moulton_koeff( int );

//Adams-Milne method
vec AdMiln_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr AdMiln_function_solve_vector ( float, float, float, vec ( * ) ( float, vec), vec );

//метотоды Адамса-Башфорта
vec adams5_function_solve ( float, float, float, int, float ( * ) ( float, float ), float);
matr adams5_function_solve_vector ( float, float, float, int, vec ( * ) ( float, vec ), vec );
matr adams5_system_solve_matrix ( float, int, vec, matr );
matr adams5_system_solve_matrix_omp ( float, int, vec, matr );
matr adams5_system_solve_matrix_mpi ( float, int, vec, matr );

//Adams-Moulton method
vec AdMltn_function_solve ( float, float, float, int, float ( * ) ( float, float ), float);
matr AdMltn_function_solve_vector ( float, float, float, int, vec ( * ) ( float, vec), vec );
matr AdMltn_system_solve_matrix ( float, int, vec, matr );
matr AdMltn_system_solve_matrix_omp ( float, int, vec, matr );
matr AdMltn_system_solve_matrix_mpi ( float, int, vec, matr );

//Euler method
vec eiler_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr eiler_system_solve_vector ( float, float, float, int, vec ( * ) ( float, vec ), vec );
matr eiler_system_solve_matrix ( float, int, vec, matr);
matr eiler_system_solve_matrix_omp ( float, int, vec, matr);
matr eiler_system_solve_matrix_mpi ( float, int, vec, matr);

//Runge-Kutta classical 4 stage method
vec rk4_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr rk4_system_function_solve_vector ( float, float, float, int, vec( * ) ( float, vec), vec );
matr rk4_system_solve_matrix ( float, int, vec, matr);
matr rk4_system_solve_matrix_omp ( float, int, vec, matr);
matr rk4_system_solve_matrix_mpi ( float, int, vec, matr);

//Runge-Kutta classical 6 stage method
vec rk6_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr rk6_system_function_solve_vector ( float, float, float, int, vec( * ) ( float, vec), vec );
matr rk6_system_solve_matrix ( float, int, vec, matr);
matr rk6_system_solve_matrix_omp ( float, int, vec, matr);
matr rk6_system_solve_matrix_mpi ( float, int, vec, matr);

//Runge-Kutta Dormand-Prince method
matr rk54( float, float, float ( * ) ( float, float ), float , float, float, float);
matr rk54( float, float, float, vec( * ) ( float, vec), vec, float , float, float );

//Runge-Kutta Fahlberg method
matr rk45( float, float, float ( * ) ( float, float ), float , float, float, float);
matr rk45( float, float, float, vec( * ) ( float, vec), vec, float , float, float );
}

#endif

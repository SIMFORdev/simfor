#include <simfor/odu.hpp>
#include <simfor/types.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <iostream>
#include <fstream>

double odu ( double t, double x )
    {
    //return x*x - x*x*x;
    //return t*t - 2*x;
    return ( 6*x - 13*t*t*t - 22*t*t + 17*t - 11 + sin ( t ) );
    }
double solution ( double t )
    {
    //return ( 13.*t*t*t ) /6 + ( 19.*t*t ) /4 - ( 5.*t ) /4 + 13./8 - cos ( t ) /37. - ( 6*sin ( t ) ) /37. + ( 119.*exp ( 6*t ) ) /296;
    return ( 3./4 * exp ( -2*t ) + 1./2*t*t - 1./2*t + 1./4 );
    }

void write_matr ( const std::string file, const simfor::matr A )
    {
    std::ofstream data ( file, std::ofstream::out );
    for ( int i = 0; i < A.size1(); ++i )
        {
        for ( int j = 0; j < A.size2(); ++j )
            data << A ( i, j ) << " ";
        data << std::endl;
        }
    }
simfor::vec absvec(simfor::vec v)
    {
    simfor::vec vabs (v.size());
    for (unsigned i = 0; i < v.size(); ++i)
        vabs(i) = std::abs ( v(i) );
    return vabs;
    }



int main()
    {

    double a = 0, b = 5, x0 = solution(a);
    int n = 100;
    double h = 1e-2;
    simfor::matr adaptive;
    simfor::matr data ( 4, n+1 );

    for ( int i = 0; i < n+1; ++i )
        {
        data ( 0, i ) = ( a + ( b-a ) / n * i );
        data ( 1, i ) = solution ( data ( 0, i ) );
        }

    row ( data, 2 ) = simfor::eiler_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    write_matr ( "odu_1_Eul.txt", data);
//     std::cout << row ( data, 2 ) ( n ) - solution ( b ) << "\n";

    row ( data, 2 ) = simfor::rk4_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    write_matr ( "odu_1_RK4.txt", data);
//     std::cout << row ( data, 2 ) ( n ) - solution ( b ) << "\n";
    adaptive = simfor::rk54 ( a, b, odu, x0, 10e-3, 10e-4, 10e-1 );
    write_matr ( "odu_1_RK54.txt", adaptive);
    std::cout << solution(b) - adaptive ( 1, adaptive.size2()-1 )<< "\n";

    row ( data, 2 ) = simfor::AdMiln_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    write_matr ( "odu_1_AdMiln.txt", data);
//     std::cout << row ( data, 2 ) ( n ) - solution ( b ) << "\n";

    row ( data, 2 ) = simfor::adams5_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    write_matr ( "odu_1_AdBash.txt", data);
//     std::cout << row ( data, 2 ) ( n ) - solution ( b ) << "\n";

    row ( data, 2 ) = simfor::AdMltn_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    write_matr ( "odu_1_AdMltn.txt", data);
//     std::cout << y ( n ) - solution ( b ) << "\n";

    return 0;
    }

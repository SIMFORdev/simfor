#include <simfor/odu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>


int main ( int argc, char* argv [] )
    {
    double a = 0, b = 5, h;
    h = 1e-5;
    int n = ceil ( (b-a)/h );
    simfor::matr F(4, 4) , Y;

    // double matr[] = {119.46, 185.38, 126.88, 121.03,-10.395,-10.136,-3.636,8.577,-53.302,-85.932,-63.182,-54.211,-115.58,-181.57,-112.8,-199};

    std::ifstream file("/home/majong/Documents/simfor/examples/BOOST_ODE/matrix.txt");
    file >> F;

    simfor::vec x0 ( 4 ) ;

    x0 ( 0 ) = x0 ( 1 ) = x0 ( 2 ) = x0 ( 3 ) = 1;



    omp_set_num_teams(8);
    double t, avg_t = 0, min_t = 1000;

    t = clock();
    //Y = simfor::eiler_system_solve_matrix_omp ( h, n, x0, F );
    t = ( clock() - t ) / CLOCKS_PER_SEC ;

    std::cout << "ESEQ " << t  << "\n";

    t = clock();
    Y = simfor::rk4_system_solve_matrix_omp ( h, n, x0, F );
    t = ( clock() - t ) / CLOCKS_PER_SEC ;
    std::cout << Y << "\n";
    std::cout << "RKSEQ " << t  << "\n";



    t = clock();
    //Y = simfor::adams5_system_solve_matrix_omp ( h, n, x0, F );
    t = ( clock() - t ) / CLOCKS_PER_SEC ;
    std::cout << "ADBSEQ " << t << "\n";

     t = clock();
   //Y = simfor::AdMltn_system_solve_matrix ( h, n, x0, F );

    t = ( clock() - t ) / CLOCKS_PER_SEC ;
    std::cout << "AdMltnSEQ " << t << "\n";

    }

#include <simfor/elementary.hpp>

namespace simfor
{
double v_norm2 ( vec a, vec b )
    {
    double tb = 0, ta = 0;

    for ( int i = 0; i < a.size (); i++ )
        {
        ta += a[i]*a[i];
        tb += b[i]*b[i];
        }
    return sqrt ( tb-ta );
    }

double v_norm2_omp ( vec a, vec b )
    {
    double tb = 0, ta = 0;
    #pragma omp parallel for private(tb, ta)
    for ( int i = 0; i < a.size (); i++ )
        {
        ta += a[i]*a[i];
        tb += b[i]*b[i];
        }
    return sqrt ( tb-ta );
    }

double m_norm4 ( matr A )
    {
    double res = 0;
    for ( int i = 0; i < A.size1(); i++ )
        for ( int j = 0; j < A.size2(); j++ )
            res += A ( j, i ) * A ( j, i );
    return sqrt ( res );
    }

double m_norm4_omp ( matr A )
    {
    double res = 0;
    #pragma omp parallel for reduction(+:res)
    for ( int i = 0; i < A.size1(); i++ )
        for ( int j = 0; j < A.size2(); j++ )
            res += A ( j, i ) *A ( j, i );
    return sqrt ( res );
    }

double psum ( vec a )
    {
    double res = 0;
    for ( int i = 0; i < a.size(); i++ )
        res += a ( i );
    return res;
    }

double psum_omp ( vec a )
    {
    double res = 0;
    double time;
    time = omp_get_wtime();
    #pragma omp parallel for reduction(+:res)
    for ( int i = 0; i < a.size(); i++ )
        res += a ( i );
    time = time - omp_get_wtime();
    return res;
    }
}

#include <iostream>
#include <cmath>
#include "simfor/ZeidelMpi.hpp"

using namespace std;

namespace simfor{

    /**
     * Parallel solution of the linear system Ax = b using the Zeidel method.
     *
     * @param mat The coefficient matrix A of the system.
     * @param vecB The right-hand side vector b of the system.
     * @param N The number of rows in the matrix A.
     * @return The solution vector x of the system Ax = b.
     */
    vec ZeidelMpi(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N);
        long double eps = 1e-6;
        size_t iter{}, iterLim = 100*N;
        long double sum = {};

        mpi::communicator world;
        int p = world.size();
        int r = world.rank();
        do{
            // Calculate x_n[i] for each processor
            int cnt = (N / p);
            int from = r * cnt;
            int to = N;
            if ( r != p-1 ) to = from + cnt;

            for(auto i = from; i < to; i++){
                x_n(i) = vecB(i);
                for(auto j = 0; j < i; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n(i) -= mat(i,j) * x_n(j);
                    };
                }

                for(auto j = i; j < N; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n(i) -= mat(i,j) * x(j);
                    };
                }

                x_n[i] /= mat(i,i);
            }

            // Send data to the root
            world.barrier();
            if (r != 0){
                world.send( 0, 1, from );
                world.send( 0, 2, to );
                world.send(0, 3, &x_n(from), ( to - from));
            }else{
                // Receive data from other processors
                for ( int i = 1; i < p; ++i ){
                    world.recv ( i, 1, from);
                    world.recv ( i, 2, to);
                    world.recv ( i, 3, &x_n(from), ( to - from));
                }
            }
            // Broadcast the result to all processors
            broadcast(world, &x_n(0), N,  0);

            // Check convergence
            from = r * cnt;
            to = N-1;
            if ( r != p-1 ) to = from + cnt;
            int tmp{0}, flag{0};
			for(int i = from; i < to; ++i){
				tmp += (fabsf64(x_n(i)-x(i)) > 1e-9 ? 1 : 0);
			}
			// Reduce the number of non-converged elements
			reduce(world, tmp, flag, std::plus<int>(), 0);
			// Broadcast the convergence flag to all processors
			broadcast(world, flag, 0);
            world.barrier();

            // Assign the new values to x
            x.assign(x_n);

            if (flag==0) break;

        } while (1);
        
        return x;
    }

}
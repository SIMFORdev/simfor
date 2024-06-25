#include <simfor/SimpleIterMpi.hpp>

namespace simfor{

    /**
     * Simple iterative method for solving linear systems with MPI parallelization
     * @param mat Coefficient matrix
     * @param vecB Right hand side vector
     * @param N Number of unknowns
     * @return Solution vector
     */
    vec SimpleIterMpi(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N); // current and previous solutions
        long double eps = 1e-9;
        size_t iter{}; // number of iterations
        long double sum{}; // sum of differences between x and x_n vectors

        mpi::communicator world;
        int p = world.size(); // number of MPI processes
        int r = world.rank(); // rank of the current process
        int cnt = (N / p); // number of unknowns per process
        int from = r * cnt; // first unknown for the current process
        int to = N; // last unknown for the current process
        if ( r != p-1 ) to = from + cnt; // if not the last process

        for (auto i = from; i < to; i++){ // initialize the solution vector
            x(i) = vecB(i) / mat(i,i);
        }

        do{ // do until convergence
            from = r * cnt; // first unknown for the current process
            to = N; // last unknown for the current process
            if ( r != p-1 ) to = from + cnt; // if not the last process

            for(auto i = from; i < to; i++){ // calculate new solution at i-th unknown
                x_n(i) = vecB(i) / mat(i,i); // initialize new solution
                for(auto j = 0; j < N; j++){ // iterate over all unknowns
                    if (i == j)
                        continue; // skip the same unknown
                    else {
                        x_n(i) -= mat(i,j) / mat(i,i) * x(j); // calculate new solution
                    };
                }
            }

            world.barrier(); // wait for all processes
            if (r != 0){ // if not the master process
                world.send( 0, 1, from ); // send first unknown index
                world.send( 0, 2, to ); // send last unknown index
                world.send(0, 3, &x_n(from), ( to - from)); // send new solution vector
            }else{ // if master process
                for ( int i = 1; i < p; ++i ){ // iterate over all processes
                    world.recv ( i, 1, from); // receive first unknown index
                    world.recv ( i, 2, to); // receive last unknown index
                    world.recv ( i, 3, &x_n(from), ( to - from)); // receive new solution vector
                }
            }
            broadcast(world, &x_n(0), N,  0); // broadcast new solution vector to all processes

            from = r * cnt; // first unknown for the current process
            to = N-1; // last unknown for the current process
            if ( r != p-1 ) to = from + cnt; // if not the last process
            int tmp{0}, flag{0}; // number of different elements and convergence flag
			for(int i = from; i < to; ++i){ // calculate number of different elements
				tmp += (fabsf64(x_n(i)-x(i)) > 1e-9 ? 1 : 0);
			}
			reduce(world, tmp, flag, std::plus<int>(), 0); // reduce number of different elements
			broadcast(world, flag, 0); // broadcast convergence flag to all processes
            world.barrier(); // wait for all processes

            x.assign(x_n); // assign new solution vector to the previous one

            if (flag==0) break; // if solution is converged, break the loop

        } while (1); // continue iterating until convergence

        return x; // return the converged solution vector
    }



}
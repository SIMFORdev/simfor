#include <simfor/TridiagonalMpi.hpp>

namespace simfor {
    /**
     * Solves the tridiagonal system Ax=b with MPI parallelization
     *
     * @param A The coefficient matrix
     * @param b The right hand side
     *
     * @return The solution vector x
     */
    vec TridiagonalMpi(matr &A, vec &b)
    {
        size_t N(A.size1()); // number of rows
        vec result(N);      // Solution vector
        vec U(N), V(N);     // Temporary vectors

        V(0) = A(0, 1) / (-A(0, 0)); // Initialization
        U(0) = -b(0) / (-A(0, 0));

        mpi::communicator world;

        int p = world.size(); // Number of processes
        int r = world.rank(); // Process rank

        int cnt = (N / p); // Number of elements per process
        int from = r * cnt; // The first index of the current process
        int to = N-1; // The last index of the current process
        if ( r != p-1 ) to = from + cnt; // If it's not the last process
        if ( r == 0) from = 1; // If it's the first process

        for (size_t i = from; i < to; i++){
            V(i) = A(i, i+1) / (-A(i,i) - A(i,i-1) * V(i-1));
            // Calculate V[i] using the result of V[i-1]

            U(i) = (A(i, i-1)*U(i-1) - b(i)) / (-A(i, i) - A(i,i-1)*V(i-1));
            // Calculate U[i] using the result of U[i-1] and V[i-1]
        }

        world.barrier(); // Wait for all processes to finish calculation
        if (r != 0){
            world.send( 0, 1, from ); // Send from to rank 0
            world.send( 0, 2, to ); // Send to to rank 0
            world.send(0, 3, &V(from), ( to - from)); // Send V[from] to rank 0
            world.send(0, 4, &U(from), ( to - from)); // Send U[from] to rank 0
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from); // Receive from process i
                world.recv ( i, 2, to); // Receive to process i
                world.recv ( i, 3, &V(from), ( to - from)); // Receive V[from] from process i
                world.recv ( i, 4, &U(from), ( to - from)); // Receive U[from] from process i
            }
        }

        broadcast(world, &V(0), N,  0); // Broadcast V to all processes
	    broadcast(world, &U(0), N,  0); // Broadcast U to all processes
        world.barrier(); // Wait for all processes to finish broadcasting

        V(N-1) = 0; // Last element of V
        U(N-1) = (A(N-1,N-2) * U(N-2) - b(N-1)) / (-A(N-1, N-1) - A(N-1, N-2) * V(N-2));
        // Calculate U[N-1]
        result(N-1) = U(N-1);
        
        for (size_t i = N-1; i > 0; i--){ //n-1 to 0
            result(i-1) = V(i-1) * result(i) + U(i-1);
        } // Calculate the solution
        broadcast(world, &result(0), N,  0); // Broadcast solution to all processes

        return result; // Return the solution vector
    }

}
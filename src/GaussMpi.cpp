#include "simfor/GaussMpi.hpp"

namespace simfor
{
    /**
     * @brief Solves linear system Ax = b using Gauss-Jordan method with MPI.
     * @param mat Coefficient matrix.
     * @param B Right-hand side vector.
     * @param res Solution vector.
     * @param N Size of the system.
     */
    void GaussianEliminationMpi(matr &mat, vec &B, vec &res, int N){
        // Local matrices and result
        matr L(N, N), U(N, N), L_inv(N, N), T(N, N);
        vec C(N);
        // MPI variables
        mpi::communicator world;
        int p = world.size(); // Number of processes
        int r = world.rank(); // Rank of the current process

        // Compute local submatrices
        int cnt = floor(N / p); // Number of rows for each process
        int from = r * cnt; // First row for current process
        int to = N; // Last row for current process
        if ( r != p-1 ) to = from + cnt; // Last row for last process

        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j)
                if (j <= i){
                    L(i, j) = mat(i, j);
                    U(i, j) = 0;
                }else{
                    U(i, j) = mat(i, j);
                    L(i, j) = 0;
                }
        }

        world.barrier(); // Synchronize processes before sending
        if (r != 0){
            world.send( 0, 1, from ); // Send first row and last row for current process
            world.send( 0, 2, to );
            world.send(0, 3, &L(from), ( to - from) * N); // Send L matrix
            world.send(0, 4, &U(from), ( to - from) * N); // Send U matrix
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from ); // Receive first row and last row from each process
                world.recv ( i, 2, to );
                world.recv ( i, 3, &L(from), ( to - from) * N); // Receive L matrix
                world.recv ( i, 4, &U(from), ( to - from) * N); // Receive U matrix
            }
        }

        // Compute inverse of L and L*U
        cnt = floor(N / p);
        from = r * cnt;
        to = N;
        if ( r != p-1 ) {to = from + cnt;}
        broadcast(world, &L.at_element(0, 0), N*N,  0); // Broadcast L matrix
        broadcast(world, &U.at_element(0, 0), N*N,  0); // Broadcast U matrix
        L_inv = findInvMatGaussJordan(L, N); // Compute L^-1

        cnt = floor(N / p);
        from = r * cnt;
        to = N;
        if ( r != p-1 ) to = from + cnt;
        matmult(L_inv, U, T, N); // Compute L^-1 * U
            
        world.barrier(); // Synchronize processes before computing C
        cnt = floor(N / p);
        from = r * cnt;
        to = N;
        if ( r != p-1 ) to = from + cnt;

        for (int i = from; i < to; ++i){
            C(i) = 0.0;
            for (int j = 0; j < N;  ++j){
                C(i) += L_inv(i, j) * B(j);
            }
        }

        if (r != 0){
            world.send(0, 1, from); // Send first row and last row
            world.send(0, 2, to);
            world.send(0, 3, &C(from), (to-from)); // Send C vector
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv(i, 1, from);
                world.recv(i, 2, to);
                world.recv(i, 3, &C(from), (to-from));
            }
        }

        world.barrier(); // Synchronize processes before iterative process
        vec resCycle(N); // Intermediate solution vector
        resCycle.clear(); // Clear NaNs
        int flag{0}; // Flag for convergence
        do{
            world.barrier(); // Synchronize processes before iteration
            flag = 0;
            res.assign(resCycle);
            cnt = floor(N / p);
            from = r * cnt;
            to = N;
            if ( r != p-1 ) {to = from + cnt;}
            for (int i = from; i < to; ++i){
                resCycle(i) = 0;
                for (int j = 0; j < N;  ++j){
                    resCycle(i) += T(i, j) * res(j);
                }
                resCycle(i) += C(i);
            }
            world.barrier(); // Synchronize processes before sending
            if (r != 0){
                world.send(0, 1, from);
                world.send(0, 2, to);
                world.send(0, 3, &resCycle(from), (to-from));
            }else{
                for (int i = 1; i < p; ++i ){
                    world.recv(i, 1, from);
                    world.recv(i, 2, to);
                    world.recv(i, 3, &resCycle(from), (to-from));
                }
            }
            broadcast(world, &resCycle(0), N,  0); // Broadcast intermediate result
            cnt = floor(N / p);
            from = r * cnt;
            to = N;
            if ( r != p-1 ) {to = from + cnt;}
            int tmp{0}; // Dummy variable for summation
            broadcast(world, &res(0), N, 0); // Broadcast final result
            for(int i = from; i < to; ++i){
                tmp += (fabsf64(resCycle(i)-res(i)) > 1e-9 ? 1 : 0);
            }
            reduce(world, tmp, flag, std::plus<int>(), 0); // Summation of convergence flag
            broadcast(world, flag, 0); // Broadcast convergence flag
            world.barrier(); // Synchronize processes
        } while (flag>1);
    }

    matr findInvMatGaussJordan(matr &mat_orig, int order){
        matr mat(order, order*2);
        matr res(order, order);
        mat.clear();	//IT CLEARS NAN NUMBERS! DO NOT DELETE!!!
        int i, j, k;
        double temp;
        mpi::communicator world;

        int p = world.size();
        int r = world.rank();

        int cnt = order / p;
        int from = r * cnt;
        int to = order;
        if ( r != p-1 ) to = from + cnt;
        world.barrier();

        for (i = from; i < to; ++i){
            for (j = 0; j < order; j++){
                mat(i,j) = mat_orig(i, j);
            }
        }

        for (i = from; i < to; i++) {
            for (j = 0; j < 2 * order; j++) {
                    if (j == (i + order)){
                        mat(i,j) = 1;
                    }
                }
        }

        world.barrier();
        if (r != 0){
            world.send(0, 1, from);
            world.send(0, 2, to);
            world.send(0, 3, &mat.at_element(from, 0), ( to - from) * order*2);
        }else{
            for (i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &mat.at_element(from, 0), ( to - from) * order*2);
            }
        }
        broadcast(world, &mat(0, 0) , order*order*2,  0);
        world.barrier();
        from = r * cnt;
        to = order;
        if ( r != p-1 ) to = from + cnt;
        
        for (i = to-1; i > from; i--) {
            if (mat(i-1,0) < mat(i,0)) {
                std::swap(mat.at_element(i, 0), mat.at_element(i - 1, 0));
            }
        }

        world.barrier();
        if (r != 0){
            world.send(0, 1, from);
            world.send(0, 2, to);
            world.send(0, 3, &mat.at_element(from, 0), ( to - from) * order*2);
        }else{
            for (i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &mat.at_element(from, 0), ( to - from) * order*2);
            }
        }

        if(r == 0) {
            for (i = 0; i < order; i++) {
                for (j = 0; j < order; j++) {
                    if (j != i) {
                        temp = mat(j,i) / mat(i,i);
                        for (k = 0; k < 2 * order; k++) {
                            mat(j,k) -= mat(i,k) * temp;
                        }
                    }
                }
            }
        }

        broadcast(world, &mat(0, 0) , order*order*2,  0);
        world.barrier();
        from = r * cnt;
        to = order;
        if ( r != p-1 ) to = from + cnt;

        for (i = from; i < to; i++) {
            temp = mat(i,i);
            for (j = 0; j < 2 * order; j++) {
                mat(i,j) = mat(i,j) / temp;
            }
        }

        if ( r != p-1 ) to = from + cnt;
        for (i = from; i < to; ++i){
                for (j = 0; j < order; j++){
                    res(i, j) = mat(i,j+order);
                }
        }
        if (r != 0){
            world.send(0, 1, from);
            world.send(0, 2, to);
            world.send(0, 3, &res.at_element(from, 0), ( to - from) * order);
        }else{
            for (i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &res.at_element(from, 0), ( to - from) * order);
            }
        }
        broadcast(world, &res(0, 0) , order*order, 0);
        return res;
    }

    void matmult( matr &A, matr &B, matr &C, int n )
        {
        mpi::communicator world;

        int p = world.size();
        int r = world.rank();

        int cnt = n / p;
        int from = r * cnt;
        int to = n;

        if ( r != p-1 ) to = from + cnt;

        for ( int i = from; i < to; ++i )
            for ( int j = 0; j < n;  ++j )
                for ( int k = 0; k < n; ++k )
                C( i, j ) += (-A( i, k )) * B( k, j );

        world.barrier ();

        if ( r != 0 )
            {
            world.send ( 0, 1, from );
            world.send ( 0, 2, to );
            world.send ( 0, 3, &C.at_element(from, 0), ( to - from ) * n );
            }
        else
            {
            for ( int i = 1; i < p; ++i )
                {
                world.recv ( i, 1, from );
                world.recv ( i, 2, to );
                world.recv ( i, 3, &C.at_element(from, 0), ( to - from ) * n );
                }
        }
    }
} // namespace simfor

#include <simfor/TridiagonalMpi.hpp>

namespace simfor {
    vec TridiagonalMpi(matr &A, vec &b){
        size_t N(A.size1());
        vec result(N);
        vec U(N), V(N);
        V(0) = A(0, 1) / (-A(0, 0));
        U(0) = -b(0) / (-A(0, 0));

        mpi::communicator world;

        int p = world.size();
        int r = world.rank();

        int cnt = (N / p);
        int from = r * cnt;
        int to = N-1;
        if ( r != p-1 ) to = from + cnt;
        if ( r == 0) from = 1;

        for (size_t i = from; i < to; i++){
            V(i) = A(i, i+1) / (-A(i,i) - A(i,i-1) * V(i-1));

            U(i) = (A(i, i-1)*U(i-1) - b(i)) / (-A(i, i) - A(i,i-1)*V(i-1));
        }
        
        world.barrier();
        if (r != 0){
            world.send( 0, 1, from );
            world.send( 0, 2, to );
            world.send(0, 3, &V(from), ( to - from));
            world.send(0, 4, &U(from), ( to - from));
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &V(from), ( to - from));
                world.recv ( i, 4, &U(from), ( to - from));
            }
        }
        broadcast(world, &V(0), N,  0);
	    broadcast(world, &U(0), N,  0);
        world.barrier();

        V(N-1) = 0;
        U(N-1) = (A(N-1,N-2) * U(N-2) - b(N-1)) / (-A(N-1, N-1) - A(N-1, N-2) * V(N-2));
        result(N-1) = U(N-1);
        
        for (size_t i = N-1; i > 0; i--){ //n-1 to 0
            result(i-1) = V(i-1) * result(i) + U(i-1);
        }
        broadcast(world, &result(0), N,  0);
        return result;
    }
}
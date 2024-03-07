#include "simfor/GaussMpi.hpp"

namespace simfor
{
    void GaussianEliminationMpi(matr &mat, vec &B, vec &res, int N){
        matr L(N, N), U(N, N), L_inv(N, N), T(N, N);
        vec C(N);
        mpi::communicator world;

        int p = world.size();
        int r = world.rank();

        int cnt = floor(N / p);
        int from = r * cnt;
        int to = N;
        if ( r != p-1 ) to = from + cnt;

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

        world.barrier();
        if (r != 0){
            world.send( 0, 1, from );
            world.send( 0, 2, to );
            world.send(0, 3, &L(from), ( to - from) * N);
            world.send(0, 4, &U(from), ( to - from) * N);
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv ( i, 1, from);
                world.recv ( i, 2, to);
                world.recv ( i, 3, &L(from), ( to - from) * N);
                world.recv ( i, 4, &U(from), ( to - from) * N);
            }
        }

        cnt = floor(N / p);
        from = r * cnt;
        to = N;
        if ( r != p-1 ) {to = from + cnt;}
        broadcast(world, &L.at_element(0, 0), N*N,  0);
        broadcast(world, &U.at_element(0, 0), N*N,  0);
        L_inv = findInvMatGaussJordan(L, N);

        cnt = floor(N / p);
        from = r * cnt;
        to = N;
        if ( r != p-1 ) to = from + cnt;
        matmult(L_inv, U, T, N);
            
        world.barrier();
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
            world.send(0, 1, from);
            world.send(0, 2, to);
            world.send(0, 3, &C(from), (to-from));
        }else{
            for ( int i = 1; i < p; ++i ){
                world.recv(i, 1, from);
                world.recv(i, 2, to);
                world.recv(i, 3, &C(from), (to-from));
            }
        }

        world.barrier();
        vec resCycle(N);
        resCycle.clear(); //IT CLEARS NAN NUMBERS! DO NOT DELETE!!!
        int flag{0};
        do{
            world.barrier();
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
            world.barrier();
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
            broadcast(world, &resCycle(0), N,  0);
            cnt = floor(N / p);
            from = r * cnt;
            to = N;
            if ( r != p-1 ) {to = from + cnt;}
            int tmp{0};
            broadcast(world, &res(0), N, 0);
            for(int i = from; i < to; ++i){
                tmp += (fabsf64(resCycle(i)-res(i)) > 1e-9 ? 1 : 0);
            }
            reduce(world, tmp, flag, std::plus<int>(), 0);
            broadcast(world, flag, 0);
            world.barrier();
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

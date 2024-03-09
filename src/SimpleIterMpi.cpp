#include "simfor/SimpleIterMpi.hpp"

namespace simfor{

    vec SimpleIterMpi(matr &mat, vec &vecB, int N){
        vec x_n(N); vec x(N);
        long double eps = 1e-9;
        size_t iter{};
        long double sum{};

        mpi::communicator world;
        int p = world.size();
        int r = world.rank();
        int cnt = (N / p);
        int from = r * cnt;
        int to = N;
        if ( r != p-1 ) to = from + cnt;

        for (auto i = from; i < to; i++){
            x(i) = vecB(i) / mat(i,i);
        }

        do{
            from = r * cnt;
            to = N;
            if ( r != p-1 ) to = from + cnt;

            for(auto i = from; i < to; i++){
                x_n(i) = vecB(i) / mat(i,i);
                for(auto j = 0; j < N; j++){
                    if (i == j)
                        continue;
                    else {
                        x_n(i) -= mat(i,j) / mat(i,i) * x(j);
                    };
                }
            }

            world.barrier();
            if (r != 0){
                world.send( 0, 1, from );
                world.send( 0, 2, to );
                world.send(0, 3, &x_n(from), ( to - from));
            }else{
                for ( int i = 1; i < p; ++i ){
                    world.recv ( i, 1, from);
                    world.recv ( i, 2, to);
                    world.recv ( i, 3, &x_n(from), ( to - from));
                }
            }
            broadcast(world, &x_n(0), N,  0);

            from = r * cnt;
            to = N-1;
            if ( r != p-1 ) to = from + cnt;
            int tmp{0}, flag{0};
			for(int i = from; i < to; ++i){
				tmp += (fabsf64(x_n(i)-x(i)) > 1e-9 ? 1 : 0);
			}
			reduce(world, tmp, flag, std::plus<int>(), 0);
			broadcast(world, flag, 0);
            world.barrier();

            x.assign(x_n);

            if (flag==0) break;

        } while (1);

        return x;
    }


}
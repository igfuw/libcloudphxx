#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then $cmake ../; fi 
$cmake -DCMAKE_BUILD_TYPE=Debug ../

#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++ -show
#ldd //usr/lib/libmpi_cxx.so
#ldd /usr/local/clang-7.0.0/lib/libomp.so
#ldd /usr/lib/x86_64-linux-gnu/libboost_mpi.so
#ldd /home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/lib/libmpi.so


## manual build of libcloud
#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++  -DUSE_MPI -Dcloudphxx_lgrngn_EXPORTS -I/usr/local/include -I/home/travis/build/igfuw/libcloudphxx/include  -fPIC   -Wextra -g -Og -DTHRUST_DEBUG -fopenmp=libomp -std=gnu++11 -o CMakeFiles/cloudphxx_lgrngn.dir/src/lib.cpp.o -c /home/travis/build/igfuw/libcloudphxx/src/lib.cpp
#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++  -DUSE_MPI -Dcloudphxx_lgrngn_EXPORTS -I/usr/local/include -I/home/travis/build/igfuw/libcloudphxx/include  -fPIC   -Wextra -g -Og -DTHRUST_DEBUG -fopenmp=libomp -std=gnu++11 -o CMakeFiles/cloudphxx_lgrngn.dir/src/lib_cpp.cpp.o -c /home/travis/build/igfuw/libcloudphxx/src/lib_cpp.cpp
#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++  -DUSE_MPI -Dcloudphxx_lgrngn_EXPORTS -I/usr/local/include -I/home/travis/build/igfuw/libcloudphxx/include  -fPIC   -Wextra -g -Og -DTHRUST_DEBUG -fopenmp=libomp -std=gnu++11 -o CMakeFiles/cloudphxx_lgrngn.dir/src/lib_omp.cpp.o -c /home/travis/build/igfuw/libcloudphxx/src/lib_omp.cpp
## link
#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++ -fPIC   -shared -Wl,-soname,libcloudphxx_lgrngn_dbg.so -o libcloudphxx_lgrngn_dbg.so CMakeFiles/cloudphxx_lgrngn.dir/src/lib.cpp.o CMakeFiles/cloudphxx_lgrngn.dir/src/lib_cpp.cpp.o CMakeFiles/cloudphxx_lgrngn.dir/src/lib_omp.cpp.o /usr/local/clang-7.0.0/lib/libomp.so -lpthread
#
#
## manual build of tests particles
#cd /home/travis/build/igfuw/libcloudphxx/build/tests/particles && /home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++   -I/usr/local/include -I/home/travis/build/igfuw/libcloudphxx/include  -fopenmp=libomp -std=gnu++11 -o CMakeFiles/test_particles.dir/tests_particles.cpp.o -c /home/travis/build/igfuw/libcloudphxx/tests/particles/tests_particles.cpp
## link
## from cmake
##/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++     CMakeFiles/test_particles.dir/tests_particles.cpp.o  -o test_particles -Wl,-rpath,/home/travis/build/igfuw/libcloudphxx/build ../../libcloudphxx_lgrngn_dbg.so /usr/local/clang-7.0.0/lib/libomp.so -lpthread /usr/lib/x86_64-linux-gnu/libboost_mpi.so /usr/lib/x86_64-linux-gnu/libboost_serialization.so
## no boost lib linking
#/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++     CMakeFiles/test_particles.dir/tests_particles.cpp.o  -o test_particles -Wl,-rpath,/home/travis/build/igfuw/libcloudphxx/build ../../libcloudphxx_lgrngn_dbg.so /usr/local/clang-7.0.0/lib/libomp.so -lpthread


VERBOSE=1 make

ldd libcloudphxx_lgrngn_dbg.so
ldd bindings/python/libcloudphxx.so
nm -gDC /home/travis/build/igfuw/libcloudphxx/deps/boost/lib/libboost_python3.so.1.65.1 
nm -gD /home/travis/build/igfuw/libcloudphxx/deps/boost/lib/libboost_python3.so.1.65.1 

OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
# make with RelWithDebInfo to have high optimization with asserts on
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=~/usr/local/ ../ 
VERBOSE=1 make 
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
make install
cd ../..
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

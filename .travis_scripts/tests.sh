#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC} -DBoost_NO_BOOST_CMAKE=ON; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../

/home/travis/build/igfuw/libcloudphxx/deps/mvapich2-2.3b/bin/mpic++ -show
ldd //usr/lib/libmpi_cxx.so
ldd /usr/local/clang-7.0.0/lib/libomp.so
ldd /usr/lib/x86_64-linux-gnu/libboost_mpi.so

VERBOSE=1 make
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=~/usr/local/ ../ 
VERBOSE=1 make 
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
make install
cd ../..
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

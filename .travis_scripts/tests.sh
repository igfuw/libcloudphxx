#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC} -DBoost_NO_BOOST_CMAKE=ON; fi
$cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
# make with RelWithDebInfo to have high optimization with asserts on
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
sudo make install
cd ../..
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

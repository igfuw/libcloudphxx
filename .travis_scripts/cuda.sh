#!/usr/bin/env sh
set -ex
# libcloudph++ 
mkdir build 
cd build
#if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=$PY_LIB -DPYTHON_INCLUDE_DIR=$PY_INC; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
cmake -DCMAKE_BUILD_TYPE=Release ../ 
VERBOSE=1 make 
cd ../..
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

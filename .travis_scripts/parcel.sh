#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
#if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC}; fi

# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make
sudo make install
cd ../..

# parcel tests for Release mode of libcloudph++
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install python-matplotlib python-gnuplot gnuplot-nox; fi
git clone --depth=1 git://github.com/igfuw/parcel.git
cd parcel
mkdir plots/outputs
py.test -v long_test
py.test -v unit_test
cd ..

# make libcloudph++ in Debug mode
sudo rm -rf libcloudphxx/build/*
cd libcloudphxx/build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
sudo make install
cd ../../

# parcel tests for Debug mode of libcloudph++
cd parcel
py.test -v unit_test_debug

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

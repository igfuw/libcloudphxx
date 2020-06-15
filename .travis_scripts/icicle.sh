#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# make with RelWithDebInfo to have high optimization with asserts on
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

git clone --depth=1 git://github.com/igfuw/libmpdataxx.git
cd libmpdataxx/libmpdata++
mkdir build
cd build
$cmake ..
sudo make install
cd ../../..

## icicle
cd libcloudphxx/models/kinematic_2D
mkdir -p build 
cd build
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make
ctest -VV -R travis # compare icicle results against reference data (done for full simulation for bulk schemes and a couple of steps for lagrangian)
cat Testing/Temporary/LastTest.log
cd ../../../..                                       
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

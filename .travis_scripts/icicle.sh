#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

# libmpdata (needed by icicle, skipping tests)
pwd
echo $TRAVIS_BUILD_DIR
. $TRAVIS_BUILD_DIR/.travis_scripts/get_libmpdata_dependencies.sh
sudo $apt_get_install hdf5-tools

git clone --depth=1 git://github.com/igfuw/libmpdataxx.git
cd libmpdataxx/libmpdata++
mkdir build
cd build
cmake ..
sudo make install
cd ../../..

## icicle
cd libcloudphxx/models/kinematic_2D
mkdir -p build 
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
if [[ $COMPILER == 'clang++' ]]; then make; fi # disable compilation on the CUDA machine with g++ - it has not enough memory to compile icicle
if [[ $COMPILER == 'clang++' ]]; then ctest -VV -R travis; fi # compare icicle results against reference data (done for full simulation for bulk schemes and a couple of steps for lagrangian)
if [[ $COMPILER == 'clang++' ]]; then cat Testing/Temporary/LastTest.log; fi
cd ../../../..                                       
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

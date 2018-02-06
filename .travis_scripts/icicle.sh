#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
#if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
# find python paths, taken from 
# https://github.com/breannansmith/scisim/blob/master/.travis.yml
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_INC=`python-config --includes | grep -o '\-I[^ ]*' | head -n 1 | cut -c 3-` ; fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_LIB=`python-config --ldflags | grep -o '\-L[^ ]*' | head -n 1 | cut -c 3- | xargs -I % find % -name libpython*.dylib` ; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC}; fi
# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

# libmpdata (needed by icicle, skipping tests)
pwd
echo $TRAVIS_BUILD_DIR
. $TRAVIS_BUILD_DIR/.travis_scripts/get_libmpdata_dependencies.sh
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libhdf5-7; fi
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install  -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0 libhdf5-dev; fi
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libgnuplot-iostream-dev libhdf5-serial-dev hdf5-tools cmake; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew tap homebrew/science; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install hdf5 --with-cxx; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then git clone --depth=1 https://github.com/dstahlke/gnuplot-iostream.git; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then sudo ln -s `pwd`/gnuplot-iostream/gnuplot-iostream.h /usr/local/include/gnuplot-iostream.h; fi

git clone --depth=1 git://github.com/igfuw/libmpdataxx.git
cd libmpdataxx/libmpdata++
mkdir build
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ..; fi
cmake ..
sudo make install
cd ../../..

## icicle
cd libcloudphxx/models/kinematic_2D
mkdir -p build 
cd build
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DBOOST_ROOT=/usr/local; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi
cmake -DCMAKE_BUILD_TYPE=Release ../ 
if [[ $CXX == 'clang++' ]]; then make; fi # disable compilation on the CUDA machine with g++ - it has not enough memory to compile icicle
if [[ $CXX == 'clang++' ]]; then ctest -VV -R travis; fi # compare icicle results against reference data (done for full simulation for bulk schemes and a couple of steps for lagrangian)
if [[ $CXX == 'clang++' ]]; then cat Testing/Temporary/LastTest.log; fi
cd ../../../..                                       
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

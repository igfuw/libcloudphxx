#!/usr/bin/env sh
set -ex

# libcloudph++ 
mkdir build 
cd build
# if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
# find python paths, taken from 
# https://github.com/breannansmith/scisim/blob/master/.travis.yml
#if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_INC=`python-config --includes | grep -o '\-I[^ ]*' | head -n 1 | cut -c 3-` ; fi
#if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_LIB=`python-config --ldflags | grep -o '\-L[^ ]*' | head -n 1 | cut -c 3- | xargs -I % find % -name libpython*.dylib` ; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC}; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
sudo make install
cd ../..

## drops.py (it is written in Python so no compilation, just unit tests)
# only on linux
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libhdf5-7; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install python-h5py; fi 
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then git clone --depth=1 git://github.com/igfuw/drops.py.git; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then cd drops.py; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then mkdir build; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then cd build; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then  cmake ..; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then make test || cat Testing/Temporary/LastTest.log /; fi # "/" intentional! (just to make cat exit with an error code)
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then  cd ../..; fi
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

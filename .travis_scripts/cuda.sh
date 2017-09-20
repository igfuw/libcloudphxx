#!/usr/bin/env sh
set -e
# libcloudph++ 
mkdir build 
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake ../; fi 
# find python paths, taken from 
# https://github.com/breannansmith/scisim/blob/master/.travis.yml
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_INC=`python-config --includes | grep -o '\-I[^ ]*' | head -n 1 | cut -c 3-` ; fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_LIB=`python-config --ldflags | grep -o '\-L[^ ]*' | head -n 1 | cut -c 3- | xargs -I % find % -name libpython*.dylib` ; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC}; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
cmake -DCMAKE_BUILD_TYPE=Release ../ 
VERBOSE=1 make 
cd ../..
set +e # see https://github.com/travis-ci/travis-ci/issues/6522

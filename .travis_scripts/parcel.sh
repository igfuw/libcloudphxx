#!/usr/bin/env sh
set -e

# libcloudph++ 
mkdir build 
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi # Travis default is not the packaged one

# find python paths, taken from 
# https://github.com/breannansmith/scisim/blob/master/.travis.yml
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_INC=`python-config --includes | grep -o '\-I[^ ]*' | head -n 1 | cut -c 3-` ; fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then PY_LIB=`python-config --ldflags | grep -o '\-L[^ ]*' | head -n 1 | cut -c 3- | xargs -I % find % -name libpython*.dylib` ; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DPYTHON_LIBRARY=${PY_LIB} -DPYTHON_INCLUDE_DIR=${PY_INC}; fi

# make libcloudph++ in Release mode
cmake -DCMAKE_BUILD_TYPE=Release ../ 
VERBOSE=1 make
sudo make install
cd ../..

# parcel tests for Release mode of libcloudph++
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install python-matplotlib python-gnuplot gnuplot-nox; fi
git clone --depth=1 git://github.com/igfuw/parcel.git
cd parcel
mkdir plots/outputs
py.test -v unit_test
py.test -v long_test
cd ..

# make libcloudph++ in Debug mode
cd libcloudphxx/build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make
sudo make install
cd ../../

# parcel tests for Debug mode of libcloudph++
cd parcel
py.test -v unit_test
py.test -v unit_test_debug

set +e # see https://github.com/travis-ci/travis-ci/issues/6522

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
# make with RelWithDebInfo to have high optimization with asserts on
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

# libmpdata (needed by icicle, skipping tests)
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libhdf5-7; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install  -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0 libhdf5-dev; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libboost-thread-dev libboost-timer-dev libboost-filesystem-dev libboost-iostreams-dev libgnuplot-iostream-dev libhdf5-serial-dev hdf5-tools cmake; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew tap homebrew/science; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install hdf5 --with-cxx; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then git clone --depth=1 https://github.com/dstahlke/gnuplot-iostream.git; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then sudo ln -s `pwd`/gnuplot-iostream/gnuplot-iostream.h /usr/local/include/gnuplot-iostream.h; fi

git clone --depth=1 git://github.com/igfuw/libmpdataxx.git
cd libmpdataxx/libmpdata++
mkdir build
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ..; fi # Travis default is not the packaged one
cmake ..
sudo make install
cd ../../..

## icicle
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libboost-program-options-dev; fi
cd libcloudphxx/tests/2D_cloud
mkdir -p build 
cd build
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then cmake .. -DBOOST_ROOT=/usr/local; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi # Travis default is not the packaged one
cmake -DCMAKE_BUILD_TYPE=Release ../ 
if [[ $CXX == 'clang++' ]]; then make; fi # disable compilation on the CUDA machine with g++ - it has not enough memory to compile icicle
if [[ $CXX == 'clang++' ]]; then ctest -VV -R travis; fi # compare icicle results against reference data (done for full simulation for bulk schemes and a couple of steps for lagrangian)
if [[ $CXX == 'clang++' ]]; then cat Testing/Temporary/LastTest.log; fi
cd ../../../..                                       
set +e # see https://github.com/travis-ci/travis-ci/issues/6522

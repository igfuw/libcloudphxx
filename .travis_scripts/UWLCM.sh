#!/usr/bin/env sh
set -ex

if [ $# -ne 1 ]; then
  echo "UWLCM.sh accepts exactly one argument"
  exit 1
fi

sudo $apt_get_install hdf5-tools

# libcloudph++ 
mkdir build 
cd build
# make with RelWithDebInfo to have high optimization with asserts on
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

# libmpdata
. $TRAVIS_BUILD_DIR/.travis_scripts/deps_install/libmpdata_dependencies.sh

git clone --depth=1 git://github.com/igfuw/libmpdataxx.git
cd libmpdataxx/libmpdata++
mkdir build
cd build
$cmake ..
sudo make install
cd ../../..

## UWLCM
git clone --depth=1 git://github.com/igfuw/UWLCM.git
cd UWLCM
. .travis_scripts/$1.sh
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

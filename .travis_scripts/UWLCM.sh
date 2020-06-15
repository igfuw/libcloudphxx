#!/usr/bin/env sh
set -ex

if [ $# -ne 1 ]; then
  echo "UWLCM.sh accepts exactly one argument"
  exit 1
fi

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

## UWLCM
# TEMP: use mpi branch from pdziekan
git clone --depth=1 --branch=mpi_up_to_date git://github.com/pdziekan/UWLCM.git
#git clone --depth=1 git://github.com/igfuw/UWLCM.git
cd UWLCM
. .travis_scripts/$1.sh
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

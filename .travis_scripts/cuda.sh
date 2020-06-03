#!/usr/bin/env sh
set -ex
# libcloudph++ 
mkdir build 
cd build
$cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make -j2
$cmake -DCMAKE_BUILD_TYPE=Release ../ 
VERBOSE=1 travis_wait 20 make -j2
cd ../..
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

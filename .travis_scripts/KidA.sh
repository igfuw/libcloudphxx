#!/usr/bin/env sh
set -ex

# libcloudph++ 
git remote show
git remote show origin
git fetch
git fetch origin
git status
git branch -v -a
git branch -r
git checkout -b kida-1d origin/kida-1d
git merge master
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
VERBOSE=1 make 
sudo make install
cd ../..

#Kid-A 1D
git clone --depth=1 git://github.com/igfuw/kid-libcloud.git;
cd kid-libcloud
chmod +x ./.travis_scripts/*
. ./.travis_scripts/lwp_test.sh
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

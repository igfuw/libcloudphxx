#!/usr/bin/env sh
set -ex

# libcloudph++ with the PR + merge kida-1d branch
git ls-remote --heads origin
git config --replace-all remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all
git branch -r
git branch -v -a
git stash
travis_wait 30 git merge --verbose --progress origin/kida-1d
git stash apply

mkdir build 
cd build
# /usr prefix is a hack for python3 to find libcloudphxx.so
# CMake should find the proper installation path for python3 packages (?)
# but what it finds is /usr/local/...
$cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../  -DCMAKE_INSTALL_PREFIX=/usr
VERBOSE=1 make 
sudo make install
cd ../..

#Kid-A 1D
git clone --depth=1 git://github.com/igfuw/kid-libcloud.git;
cd kid-libcloud
chmod +x ./.travis_scripts/*
. ./.travis_scripts/lwp_test.sh
set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

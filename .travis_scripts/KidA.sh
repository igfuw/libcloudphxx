#!/usr/bin/env sh
set -ex

# libcloudph++ with the PR + merge kida-1d branch
local build_head=$(git rev-parse HEAD)
git ls-remote --heads origin
git config --replace-all remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch --all
git branch -r
git branch -v -a
git status
git checkout kida-1d
git checkout ${build_head}
git merge kida-1d
#git merge origin/kida-1d
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

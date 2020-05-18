#!/usr/bin/env sh
set -e
mkdir local_clang
cd local_clang
export CLANG_VER=clang+llvm-6.0.0-x86_64-apple-darwin
wget http://releases.llvm.org/6.0.0/$CLANG_VER.tar.xz
tar -xf $CLANG_VER.tar.xz
export LLVM_PATH=`pwd`/$CLANG_VER
export CXX=$LLVM_PATH/bin/clang++
export LDFLAGS=$LDFLAGS" -L "$LLVM_PATH/lib" -Wl,-rpath,"$LLVM_PATH/lib
export CPPFLAGS=$CPPFLAGS" -I"$LLVM_PATH/include" -I"$LLVM_PATH/include/c++/v1/
cd ..

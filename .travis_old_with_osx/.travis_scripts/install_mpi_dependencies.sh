#!/usr/bin/env sh
set -e
############################################################################
## All the cached dependencies are installed in ${TRAVIS_BUILD_DIR}/deps/
#############################################################################
DEPS_DIR="${TRAVIS_BUILD_DIR}/deps"

# get libclang-dev for headers
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' ]]; then sudo $apt_get_install libclang-dev; fi
#if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' ]]; then export CXXFLAGS="-nostdinc++ ${CXXFLAGS}"; fi

# redefine CXX to the actual version used
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' ]]; then export CXX=clang++; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'g++'     ]]; then export CXX=g++; fi
# downloads and setups local clang on osx
if [[ $TEST_SUITE == 'osx_local_clang' ]]; then . ./.travis_scripts/setup_local_clang.sh; fi

#<<<<<<< HEAD
    # add a definition -DBOOST_HAS_INT128=1 to clang calls on linux to avoid errors with boost.atomic (https://svn.boost.org/trac/boost/ticket/9610)
#    if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then mkdir /tmp/bin; fi
#    if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then printf "#!/bin/sh\nexec /usr/bin/clang++ -DBOOST_HAS_INT128=1 \"\$@\"" > /tmp/bin/clang++; fi
#    if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then chmod +x /tmp/bin/clang++; fi
#    if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then sudo ln -sf /tmp/bin/clang++ /usr/bin/clang++; fi
    # put /usr/bin first to use clang++-3.5 instead of the default 3.4
    #if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then export PATH=/usr/bin:$PATH; fi

#if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then export CXXFLAGS="-DBOOST_HAS_INT128=1 ${CXXFLAGS}"; fi

# cmake 
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then wget https://github.com/Kitware/CMake/releases/download/v3.13.2/cmake-3.13.2-Linux-x86_64.sh; fi
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo sh cmake-3.13.2-Linux-x86_64.sh --prefix=/usr/local --exclude-subdir; fi

# MPI - mvapich2-2.3b
ls -A ${DEPS_DIR}/mvapich2-2.3b
if [[ -z "$(ls -A ${DEPS_DIR}/mvapich2-2.3b)" ]]; then
  wget http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.3b.tar.gz;
  tar xf mvapich2-2.3b.tar.gz;
  cd mvapich2-2.3b;
  if [[ $COMPILER == 'g++' ]]; then ./configure --disable-fortran --enable-cxx --enable-threads=multiple --with-device=ch3:sock CC=gcc CXX=g++ --prefix=${DEPS_DIR}/mvapich2-2.3b ; fi 
  if [[ $COMPILER == 'clang++' ]]; then ./configure --disable-fortran --enable-cxx --enable-threads=multiple --with-device=ch3:sock CC=clang CXX=clang++ --prefix=${DEPS_DIR}/mvapich2-2.3b ; fi 
  make -j4;  
  make install;
  cd ..;
else
  echo "Using cached mvapich2."
fi

export PATH=${DEPS_DIR}/mvapich2-2.3b/bin:${PATH}
# LIBRARY_PATH for clang?osx?
export LD_LIBRARY_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LD_LIBRARY_PATH}
export LD_RUN_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LD_RUN_PATH}
export LIBRARY_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LIBRARY_PATH}

export CXX=${DEPS_DIR}/mvapich2-2.3b/bin/mpic++  # full path, since libtool in hdf5 installation does not understand PATH set above (?)
export CC=${DEPS_DIR}/mvapich2-2.3b/bin/mpicc 


# numpy needs to be installed before building boost python in order to build boost numpy
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then export apt_get_install="apt-get install -t xenial --no-install-recommends -y"; fi
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install python3-numpy; fi

# for MPI we need boost>=1.59 with mpi support, boost installation based on https://github.com/boostorg/compute/blob/master/.travis.yml
  if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then 
    ls -A ${DEPS_DIR}/boost
    if [[ -z "$(ls -A ${DEPS_DIR}/boost)" ]]; then
      wget http://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz 
      tar xf boost_1_65_1.tar.gz
      cd boost_1_65_1
      # configure and install
      echo "using python : 3.5 : /opt/pyenv/shims/python3 ;" >> $HOME/user-config.jam
      if [[ $COMPILER == 'g++' ]]; then echo "using gcc : 5.4 : g++ ;" > $HOME/user-config.jam; fi
      echo "using mpi : $CC ;" >> $HOME/user-config.jam
      cat $HOME/user-config.jam
      if [[ $COMPILER == 'g++' ]]; then
        ./bootstrap.sh --prefix=${DEPS_DIR}/boost/ --with-libraries=serialization,mpi,thread,date_time,system,iostreams,timer,filesystem,program_options,python
        ./b2 -d0 install
      fi
      if [[ $COMPILER == 'clang++' ]]; then 
        #clang installation taken from https://gist.github.com/jimporter/10442880
        ./bootstrap.sh --prefix=${DEPS_DIR}/boost/ --with-libraries=serialization,mpi,thread,date_time,system,iostreams,timer,filesystem,program_options,python --with-toolset=clang
        ./b2 clean
        ./b2 toolset=clang cxxflags="-std=c++14 -stdlib=libc++" linkflags="-stdlib=libc++" --prefix=${DEPS_DIR}/boost/ -j 4 stage release
        ./b2 install toolset=clang cxxflags="-std=c++14 -stdlib=libc++" linkflags="-stdlib=libc++" --prefix=${DEPS_DIR}/boost/
      fi
      cd ..
    else
      echo "Using cached boost."
    fi
    export BOOST_ROOT=${DEPS_DIR}/boost
    export LD_LIBRARY_PATH=${DEPS_DIR}/boost/lib:${LD_LIBRARY_PATH}
    export LD_RUN_PATH=${DEPS_DIR}/boost/lib:${LD_RUN_PATH}
    export LIBRARY_PATH=${DEPS_DIR}/boost/lib:${LIBRARY_PATH}
    export CPATH=${DEPS_DIR}/boost/include:${CPATH}
  fi

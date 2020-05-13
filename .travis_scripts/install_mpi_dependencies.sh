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
        ./b2 toolset=clang cxxflags="-std=c++14 -stdlib=libstdc++" linkflags="-stdlib=libstdc++" --prefix=${DEPS_DIR}/boost/ -j 4 stage release
        ./b2 install toolset=clang cxxflags="-std=c++14 -stdlib=libstdc++" linkflags="-stdlib=libstdc++" --prefix=${DEPS_DIR}/boost/
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

#!/usr/bin/env sh
set -e
############################################################################
## All the cached dependencies are installed in ${TRAVIS_BUILD_DIR}/deps/
#############################################################################
DEPS_DIR="${TRAVIS_BUILD_DIR}/deps"

. ${PWD}/.travis_scripts/install_mpi_dependencies.sh

# numpy needs to be installed before building boost python in order to build boost numpy
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then export apt_get_install="apt-get install -t xenial --no-install-recommends -y"; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install python-numpy; fi

# for MPI we need boost>=1.59 with mpi support, boost installation based on https://github.com/boostorg/compute/blob/master/.travis.yml
  if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then 
    ls -A ${DEPS_DIR}/boost
    if [[ -z "$(ls -A ${DEPS_DIR}/boost)" ]]; then
      wget http://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz 
      tar xf boost_1_65_1.tar.gz
      cd boost_1_65_1
      # configure and install
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

# Ubuntu dependency issue fix
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0; fi

# C++ support missing in Debian package ...
#if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then sudo $apt_get_install libhdf5-openmpi-dev; fi 
# ... so we are installing it manually:
  if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then 
    ls -A ${DEPS_DIR}/hdf5
    if [[ -z "$(ls -A ${DEPS_DIR}/hdf5)" ]]; then
      wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar
      tar xf hdf5-1.10.5.tar
      cd hdf5-1.10.5
      CXXFLAGS=-w CFLAGS=-w ./configure --enable-parallel --enable-cxx --enable-unsupported --enable-threadsafe --prefix=${DEPS_DIR}/hdf5/
      make
      sudo make install
      cd ..
    else
      echo "Using cached hdf5."
    fi
    export HDF5_ROOT=${DEPS_DIR}/hdf5
    export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}
    export LD_RUN_PATH=${HDF5_ROOT}/lib:${LD_RUN_PATH}
    export LIBRARY_PATH=${HDF5_ROOT}/lib:${LIBRARY_PATH}
    export CPATH=${HDF5_ROOT}/include:${CPATH}
    export PATH=${HDF5_ROOT}/bin:${PATH}
  fi

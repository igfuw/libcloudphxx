#!/usr/bin/env sh
set -e
############################################################################
## All the cached dependencies are installed in ${TRAVIS_BUILD_DIR}/deps/
#############################################################################
DEPS_DIR="${TRAVIS_BUILD_DIR}/deps"

# Ubuntu dependency issue fix
sudo $apt_get_install -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0

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

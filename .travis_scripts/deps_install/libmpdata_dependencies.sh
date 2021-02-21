#!/usr/bin/env sh
set -e

# Ubuntu dependency issue fix
sudo $apt_get_install  -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0

sudo $apt_get_install gnuplot-nox

# gnuplot-iostream
git clone -n https://github.com/dstahlke/gnuplot-iostream && cd gnuplot-iostream
git checkout 8b6e30c6ea5ee4f07ccf90d858d35a10cf67a3e2 # later commits require c++17
sudo cp gnuplot-iostream.h /usr/local/include/gnuplot-iostream.h
cd ..

# hdf5
if [[ $MPI == 'none' ]]; then sudo $apt_get_install libhdf5-dev; fi

__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $MPI != 'none' ]]; then source ${__dir}/mpi_hdf5.sh ; fi

sudo $apt_get_install hdf5-tools

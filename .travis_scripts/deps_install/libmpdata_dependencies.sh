#!/usr/bin/env sh
set -e

sudo $apt_get_install gnuplot-nox

# gnuplot-iostream
git clone -n https://github.com/dstahlke/gnuplot-iostream && cd gnuplot-iostream
git checkout 8b6e30c6ea5ee4f07ccf90d858d35a10cf67a3e2 # later commits require c++17
sudo cp gnuplot-iostream.h /usr/local/include/gnuplot-iostream.h
cd ..

sudo $apt_get_install  -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0 libhdf5-dev

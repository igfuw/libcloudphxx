#!/usr/bin/env sh
set -e

# gnuplot-iostream
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install gnuplot-nox; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install gnuplot; fi
sudo wget -O /usr/local/include/gnuplot-iostream.h https://raw.githubusercontent.com/dstahlke/gnuplot-iostream/master/gnuplot-iostream.h

#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libhdf5-7; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install  -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0 libhdf5-dev; fi
#if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libhdf5-serial-dev hdf5-tools cmake; fi #libgnuplot-iostream-dev 
# https://github.com/travis-ci/travis-ci/issues/8826
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install hdf5; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then git clone --depth=1 https://github.com/dstahlke/gnuplot-iostream.git; fi
#if [[ $TRAVIS_OS_NAME == 'osx' ]]; then sudo ln -s `pwd`/gnuplot-iostream/gnuplot-iostream.h /usr/local/include/gnuplot-iostream.h; fi

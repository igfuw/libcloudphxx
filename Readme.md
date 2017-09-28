libcloudph++ - a cloud (micro)physics library 
=======================================================================

To get more information on libcloudph++, please check: 
  - http://libcloudphxx.igf.fuw.edu.pl/
  - http://arxiv.org/abs/1310.1905
  - http://www.geosci-model-dev-discuss.net/7/8275/2014/

Compilation of libcloudph++ requires:
- a C++11 compliant compiler (optionally with OpenMP support)
- a CUDA compiler (optional)
- Thrust C++ library
- Boost C++ libraries

Compilation of the Python bindings for libcloudph++ (enabled
by default) further requires:
- Python interpretter and development files
- NumPy Python package
- Boost.Python C++ library
- Blitz++ C++ library

Compilation of the ``icicle'' test program described in the 
2015 GMD paper on libcloudph++ (and generation of all results
presented in the paper) further requires:
- libmpdata++ C++ library and its dependencies including:
  - HDF5 library with its optional C++ bindings
  - gnuplot-iostream C++ library and gnuplot

During development of libcloudph++, we are continuously testing
the code on Linux using GCC and LLVM/Clang as well as on OSX
using Apple/Clang - these are considered the supported platforms.

Compilation and execution of the examples shipped with libcloudph++ 
is easiest done using CMake, and the following instructions assume
you're using CMake. Some hints on CMake usage are included at the
end of this file.

The .travis.yml file shipped with the library contains a complete
set of commands needed to build and execute all tests programs
shipped with libcloudph++ on fresh Ubuntu and OSX installations -
it may contain useful information on obtaining the dependencies.

1. To check the dependencies and compile the library, please try:

  $ mkdir build

  $ cd build

  $ cmake ..

  $ make

  $ cd ..
  
The next two steps are optional test. Running the tests is highly
recommended to verify if the library works correctly in your 
environment. Nevertheless, in principle you can skip to step four
and install the library right away.
  
2. To perform unit and some other quick tests, please try:

  $ cd build/tests

  $ make test

  $ cd ../..

These tests should complete in a few minutes.

3. To reproduce all results from the GMD paper, please try:

  $ cd models/kinematic_2D

  $ mkdir build 

  $ cd build

  $ cmake ..

  $ make

  $ make test     

  $ cd ../../..

This can take over an hour if a GPU is available or longer if using
CPU only. 

4. To install the library system-wide, please try:

  $ cd build

  $ sudo make install

This will copy the libcloudph++ headers into the system include path
(e.g. /usr/include/libcloudph++), copy the libcloudph++-config.cmake 
file into the system share directory (e.g. /usr/share/libcloudph++) 
and install the shared library file (e.g. /usr/lib).

Some CMake hints:
- to point CMake to a non-default C++ compiler (e.g. clang++):
  $ cmake .. -DCMAKE_CXX_COMPILER=clang++ 

- to alter the installation prefix (e.g. /usr/ instead of /usr/local):
  $ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/usr

- to switch between debug and release (default) compilation modes 
  (has to be done after compiler choice):
  $ cmake .. -DCMAKE_BUILD_TYPE=Debug
  $ cmake .. -DCMAKE_BUILD_TYPE=Release
  
- two alternative ways of cleaning leftovers from a previous build 
  (including CMake cache files):
  $ rm -rf build/CMakeCache.txt build/CMakeFiles
  $ rm -rf build; mkdir build

- the output of commands executed by "make test" can be viewed with:
  $ less Testing/Temporary/LastTest.log


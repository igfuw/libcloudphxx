# Installation Guide

## Dependencies

-  C++11 or higher compiler
- Thrust C++ library
-  Boost C++ library
- For Python bindings: Python interpreter with NumPy, Blitz++ C++ library
- OpenMP support and CUDA compiler (optional)

## Building

```bash
git clone git@github.com:your-repo/libcloudphxx.git
cd libcloudphxx
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++
make install
```

CMake options:
- Selecting the compilation mode: -DCMAKE_BUILD_TYPE = Release / Debug / RelWithDebInfo 
- Selecting the compiler: -DCMAKE_CXX_COMPILER = g++ / clang++ / ...
- Choosing a custom directory for the installation: -DCMAKE_INSTALL_PREFIX = `/usr/local`, `/home/builds`, etc.
- Pointing to the location of a dependency, for example Thrust: -DTHRUST_INCLUDE_DIR = `/usr/local`
- Running the compilation in parallel for speedup: `make -jN install`, where N is the number of cores.

Conflicting versions of Boost and Thrust may cause the build to fail. If that happens, try e.g. Thrust 12.9 and Boost 1.83, with the following fix, or use Apptainer/Singularity instead of installing the dependencies manually.
```bash
git clone -b Thrust_2 --depth=1 https://github.com/pdziekan/odeint.git
cp -r odeint/include/boost/numeric/odeint/* /usr/include/boost/numeric/odeint/
rm -rf odeint
```

## Installation with Apptainer/Singularity

Apptainer is a container platform designed for scientific computing.
It allows packaging the software environment — including dependencies into a single, portable Singularity Image File (.sif).
Instructions for installation with Apptainer are available in [the UWLCM documentation](https://github.com/igfuw/UWLCM/blob/master/docs/Installation_guide.md).

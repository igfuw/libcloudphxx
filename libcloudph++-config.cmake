# needed for the OpenMP test to work in C++-only project
# (see http://public.kitware.com/Bug/view.php?id=11910)
cmake_minimum_required(VERSION 2.8.8) 

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
if(CMAKE_VERSION VERSION_GREATER 2.9)
  cmake_policy(VERSION 3.0)
endif()

############################################################################################
# the following variables will be set:
set(libcloudphxx_FOUND False)
set(libcloudphxx_INCLUDE_DIRS "")
set(libcloudphxx_LIBRARIES "cloudphxx_lgrngn")
set(libcloudphxx_CXX_FLAGS_DEBUG "")
set(libcloudphxx_CXX_FLAGS_RELEASE "")

############################################################################################
# debug mode compiler flags
set(libcloudphxx_CXX_FLAGS_DEBUG "${libcloudphxx_CXX_FLAGS_DEBUG} -std=c++11 -DBZ_DEBUG -DTHRUST_DEBUG -g") #TODO: -Og if compiler supports it?

############################################################################################
# release mode compiler flags
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
  set(libcloudphxx_CXX_FLAGS_RELEASE "${libcloudphxx_CXX_FLAGS_RELEASE} -std=c++11 -DNDEBUG -Ofast -march=native")
endif()

############################################################################################
# Boost libraries
find_package(Boost QUIET)
if(Boost_FOUND)
  set(libcloudphxx_LIBRARIES "${libcloudphxx_LIBRARIES};${Boost_LIBRARIES}")
  set(libcloudphxx_INCLUDE_DIRS "${libcloudphxx_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
else()
#TODO: check separately for optional and mandatory components
message(FATAL_ERROR "Boost (or some of its components) not found.

* To insall Boost, please try:
*   Debian/Ubuntu: sudo apt-get install libboost-all-dev
*   Fedora: sudo yum install boost-devel
")
endif()

############################################################################################
#list(REMOVE_DUPLICATES libcloudpxxx_INCLUDE_DIRS)
list(REMOVE_ITEM libcloudphxx_INCLUDE_DIRS "")
set(libcloudphxx_FOUND TRUE)

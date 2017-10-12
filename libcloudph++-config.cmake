# needed for the OpenMP test to work in C++-only project
# (see http://public.kitware.com/Bug/view.php?id=11910)
cmake_minimum_required(VERSION 2.8.8) 

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
# - CMP0060 - make CMake always keep absoult RPATHs, even if installing in implicit directory
if(CMAKE_VERSION VERSION_GREATER 2.9)
  cmake_policy(VERSION 3.0)
endif()

set(CMAKE_MACOSX_RPATH ON) # explicit, since policy CMP0042 didn't work...

############################################################################################
# the following variables will be set:
set(libcloudphxx_FOUND False)
set(libcloudphxx_INCLUDE_DIRS "")
set(libcloudphxx_LIBRARIES "")
set(libcloudphxx_CXX_FLAGS_DEBUG "")
set(libcloudphxx_CXX_FLAGS_RELWITHDEBINFO "")
set(libcloudphxx_CXX_FLAGS_RELEASE "")


############################################################################################
# libcloudphxx libs and headers 
# also work for non-default install location (i.e. for make DESTDIR=<dir> install)
set(libcloudphxx_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../include/")
set(libcloudphxx_LIBRARIES "${CMAKE_CURRENT_LIST_DIR}/../../lib/libcloudphxx_lgrngn.so")

############################################################################################
# debug mode compiler flags
set(libcloudphxx_CXX_FLAGS_DEBUG "${libcloudphxx_CXX_FLAGS_DEBUG} -std=c++11 -g -DTHRUST_DEBUG") #TODO: -Og if compiler supports it?

############################################################################################
# release with debug info mode compiler flags
set(libcloudphxx_CXX_FLAGS_RELWITHDEBINFO "${libcloudphxx_CXX_FLAGS_RELWITHDEBINFO} -std=c++11 -Ofast -march=native")

############################################################################################
# release mode compiler flags
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
  set(libcloudphxx_CXX_FLAGS_RELEASE "${libcloudphxx_CXX_FLAGS_RELEASE} -std=c++11 -DNDEBUG -Ofast -march=native -Winline")
endif()

############################################################################################
# Boost libraries
find_package(Boost)
if(Boost_FOUND)
#TODO: if boost is not linked in some program, link boost libs to libcloudphxx_lgrngn.so ?
#  set(libcloudphxx_LIBRARIES "${libcloudphxx_LIBRARIES};${Boost_LIBRARIES}")
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

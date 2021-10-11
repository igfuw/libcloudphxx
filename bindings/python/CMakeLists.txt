add_library(cloudphxx SHARED lib.cpp)
add_dependencies(cloudphxx git_revision.h)
set_target_properties(cloudphxx PROPERTIES SUFFIX ".so") # e.g. Mac defaults to .dylib which is not looked for by Python

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

target_include_directories(cloudphxx PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(cloudphxx ${Python3_LIBRARIES})

# informing the Python bindings where to find Boost.Python 
#if(${Boost_MINOR_VERSION} LESS 67)
#  message(STATUS "boost python as python")
#  find_package(Boost COMPONENTS python)
#else()
#  message(STATUS "boost python as python27")
#  find_package(Boost COMPONENTS python27)
#endif()

#find_package(Boost COMPONENTS python3)
#if (NOT Boost_FOUND)
#  message(FATAL_ERROR "
#    Boost.Python3 not found. 
#    Please install it (e.g. sudo apt-get install libboost-python-dev).
#  ")
#endif()
#target_link_libraries(cloudphxx ${Boost_LIBRARIES})

#find_package(Boost COMPONENTS numpy3)
#if (NOT Boost_FOUND)
#  message(FATAL_ERROR "
#    Boost.Numpy3 not found. 
#    Please install it (e.g. sudo apt-get install libboost-python-dev).
#  ")
#endif()
#target_link_libraries(cloudphxx ${Boost_LIBRARIES})
#target_compile_options(cloudphxx PRIVATE  -DBPNUMPY)

# boost python 1.65.0 replaced the booost::python::numeric API with boost::python::numpy
if(${Boost_MINOR_VERSION} LESS 65)
  message(STATUS "boost numpy as DBPNUMERIC")
  target_compile_options(cloudphxx PRIVATE -DBPNUMERIC)
  find_package(Boost COMPONENTS python3 REQUIRED)
  target_link_libraries(cloudphxx ${Boost_LIBRARIES})
else()
#  message(STATUS "boost numpy as numpy3")
#  find_package(Boost COMPONENTS numpy3 REQUIRED) # also finds python3
  message(STATUS "boost numpy as numpy")
  find_package(Boost COMPONENTS numpy REQUIRED) # also finds python3
  target_link_libraries(cloudphxx ${Boost_LIBRARIES})
  target_compile_options(cloudphxx PRIVATE  -DBPNUMPY)
#  if(${Boost_MINOR_VERSION} LESS 67)
#    message(STATUS "boost numpy as numpy")
#    find_package(Boost COMPONENTS numpy REQUIRED)
#  else()
#    #message(STATUS "boost numpy as numpy27")
#    #find_package(Boost COMPONENTS numpy27 REQUIRED)
#    message(STATUS "boost numpy as numpy")
#    find_package(Boost COMPONENTS numpy REQUIRED)
#  endif()
#  target_link_libraries(cloudphxx ${Boost_LIBRARIES})
#  target_compile_options(cloudphxx PRIVATE  -DBPNUMPY)
endif()

# requireing Blitz++
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #include <blitz/array.h>
    int main() {}
  " BLITZ_FOUND)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "
    Blitz++ library not found. 
    Please install it (e.g. sudo apt-get install libblitz0-dev).
  ")
else()
  target_link_libraries(cloudphxx blitz)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(cloudphxx PRIVATE -DBZDEBUG)
  endif()
endif()

target_link_libraries(cloudphxx cloudphxx_lgrngn)

#to retain rpath to libcloudphxx_lgrngn.so linked by libcloudphxx.so after installation
set_property(TARGET cloudphxx PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set_property(TARGET cloudphxx PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

# Python_SITELIB is an absolute path, we manually prepend install prefix
install ( TARGETS cloudphxx
   LIBRARY
     DESTINATION "${CMAKE_INSTALL_PREFIX}${Python3_SITELIB}"
     COMPONENT library
)

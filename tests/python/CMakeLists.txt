#TODO: search pytest?
#TODO: local search path for python modules (to make it work before make install)

# find_interp first, find_libs should find a matching version of libs
find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs)


if (NOT PYTHON_VERSION_STRING VERSION_EQUAL PYTHONLIBS_VERSION_STRING)
  message(STATUS "PYTHON_EXECUTABLE: " ${PYTHON_EXECUTABLE})
  message(STATUS "PYTHON_VERSION: " ${PYTHON_VERSION_STRING})
  message(STATUS "PYTHON_LIBRARIES: " ${PYTHON_LIBRARIES})
  message(STATUS "PYTHONLIBS_VERSION: " ${PYTHONLIBS_VERSION_STRING})
  message(FATAL_ERROR "Version mismatch between python interpreter and libraries")
endif()

# requireing NumPy
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy" ERROR_QUIET RESULT_VARIABLE NUMPY_FOUND)
if (NUMPY_FOUND GREATER 0)
  message(FATAL_ERROR "
    NumPy not found. 
    Please install it (e.g. sudo apt-get install python-numpy).
  ")
endif()

add_subdirectory(unit)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  add_subdirectory(physics)
endif()

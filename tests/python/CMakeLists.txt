add_library(cloudphxx SHARED lib.cpp)
set_target_properties(cloudphxx PROPERTIES SUFFIX ".so") # e.g. Mac defaults to .dylib which is not looked for by Python

#TODO: search for NumPy, Blitz, pytest

# informing the Python bindings where to find Python
find_package(PythonLibs REQUIRED)
target_include_directories(cloudphxx PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(cloudphxx ${PYTHON_LIBRARIES})

# informing the Python bindings where to find Boost.Python 
find_package(Boost COMPONENTS python REQUIRED QUIET)
target_link_libraries(cloudphxx ${Boost_LIBRARIES})

# ... where to find libclouph++
target_link_libraries(cloudphxx cloudphxx_lgrngn)

# looking for the Python interpretter
find_package(PythonInterp REQUIRED)
if (NOT PYTHON_VERSION_STRING VERSION_EQUAL PYTHONLIBS_VERSION_STRING)
  message(STATUS "PYTHON_EXECUTABLE: " ${PYTHON_EXECUTABLE})
  message(STATUS "PYTHON_VERSION: " ${PYTHON_VERSION_STRING})
  message(STATUS "PYTHON_LIBRARIES: " ${PYTHON_LIBRARIES})
  message(STATUS "PYTHONLIBS_VERSION: " ${PYTHONLIBS_VERSION_STRING})
  message(FATAL_ERROR "Version mismatch between python interpreter and libraries")
endif()

# non-pytest tests
foreach(test api_blk_1m api_blk_2m api_lgrngn api_common)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/python" 
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
  )
endforeach()

# pytest tests run with "python -m pytest"
foreach(test pytest_hello)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/python" 
    COMMAND ${PYTHON_EXECUTABLE} "-m" "pytest" "-s" "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
  )
endforeach()

# installation (see http://www.cmake.org/pipermail/cmake/2013-March/054138.html)
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.PREFIXES[-1])"
   OUTPUT_VARIABLE PYTHON_PREFIX
)
file ( TO_CMAKE_PATH "${PYTHON_PREFIX}" PYTHON_PREFIX )
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.getsitepackages()[-1])"
   OUTPUT_VARIABLE PYTHON_SITE_DIR
)
file ( TO_CMAKE_PATH "${PYTHON_SITE_DIR}" PYTHON_SITE_DIR )
string ( REGEX REPLACE "^${PYTHON_PREFIX}/" ""
   PYTHON_SITE_DIR "${PYTHON_SITE_DIR}"
)

install ( TARGETS cloudphxx
   LIBRARY
     DESTINATION ${PYTHON_SITE_DIR}
     COMPONENT library
)

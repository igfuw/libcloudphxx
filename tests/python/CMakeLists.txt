#TODO: search pytest?
#TODO: local search path for python modules (to make it work before make install)

find_package(PythonLibs)

find_package(PythonInterp REQUIRED)

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

# non-pytest tests
foreach(test api_blk_1m api_blk_2m api_lgrngn api_common unit_lgrngn_adve)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
  )
endforeach()

## pytest tests run with "python -m pytest"
#foreach(test pytest_hello)
#  #TODO: indicate that tests depend on the lib
#  add_test(
#    NAME ${test}
#    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
#    COMMAND ${PYTHON_EXECUTABLE} "-m" "pytest" "-s" "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
#  )
#endforeach()

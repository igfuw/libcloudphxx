add_library(cloudphxx SHARED lib.cpp)
set_target_properties(cloudphxx PROPERTIES SUFFIX ".so") # e.g. Mac defaults to .dylib which is not looked for by Python

# informing the Python bindings where to find Python
find_package(PythonLibs)
target_include_directories(cloudphxx PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(cloudphxx ${PYTHON_LIBRARIES})

# informing the Python bindings where to find Boost.Python 
find_package(Boost COMPONENTS python REQUIRED QUIET)
target_link_libraries(cloudphxx ${Boost_LIBRARIES})

# 
find_package(PythonInterp)

#
foreach(test api_blk_1m api_blk_2m api_lgrngn)
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/python" 
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
  )
endforeach()

# pytest tests run with "python -m pytest"
foreach(test pytest_hello blk_1m_pytest)
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/python" 
    COMMAND ${PYTHON_EXECUTABLE} "-m" "pytest" "-s" "${CMAKE_SOURCE_DIR}/tests/python/${test}.py"
  )
endforeach()

#TODO: indicate that tests depend on the lib

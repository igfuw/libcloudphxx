add_library(cloudphxx SHARED lib.cpp)

# informing the Python bindings where to find Python
find_package(PythonLibs)
target_include_directories(cloudphxx PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(cloudphxx ${PYTHON_LIBRARIES})

# informing the Python bindings where to find Boost.Python 
find_package(Boost COMPONENTS python REQUIRED QUIET)
target_link_libraries(cloudphxx ${Boost_LIBRARIES})

# 
find_package(PythonInterp)
add_test(
  NAME test
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/python" 
  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/test.py"
)
#TODO: indicate that tests depend on the lib

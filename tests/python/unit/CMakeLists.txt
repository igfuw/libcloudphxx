# non-pytest tests
foreach(test api_blk_1m api_blk_2m api_lgrngn api_common segfault_20150216 test_col_kernels test_uniform_init)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/unit/${test}.py"
  )
endforeach()

## pytest tests run with "python -m pytest"
foreach(test test_lgrngn_adve)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${PYTHON_EXECUTABLE} "-m" "pytest" "-s" "${CMAKE_SOURCE_DIR}/tests/python/unit/${test}.py"
  )
endforeach()

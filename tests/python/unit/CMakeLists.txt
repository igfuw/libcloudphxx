# non-pytest tests
foreach(test api_blk_1m api_blk_2m api_lgrngn api_common segfault_20150216 test_col_kernels test_terminal_velocities test_SD_removal test_uniform_init test_source test_chem_coal test_sstp_cond test_adve_scheme)
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

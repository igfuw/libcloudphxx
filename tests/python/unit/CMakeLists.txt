# non-pytest tests
foreach(test api_blk_1m api_blk_2m api_lgrngn api_common segfault_20150216 col_kernels terminal_velocities uniform_init source sstp_cond multiple_kappas adve_scheme lgrngn_subsidence sat_adj_blk_1m diag_incloud_time relax)

  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/unit/${test}.py"
  )
endforeach()

# tests with chemistry, chemistry does not work with MPI
if(NOT USE_MPI)
  foreach(test SD_removal chem_coal)
    #TODO: indicate that tests depend on the lib
    add_test(
      NAME ${test}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
      COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/unit/${test}.py"
    )
  endforeach()
endif()

## pytest tests run with "python -m pytest"
foreach(test lgrngn_adve)
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${Python3_EXECUTABLE} "-m" "pytest" "-s" "${CMAKE_SOURCE_DIR}/tests/python/unit/${test}.py"
  )
endforeach()

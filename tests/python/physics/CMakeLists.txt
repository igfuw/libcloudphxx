#TODO: search scipy?

# non-pytest tests
foreach(test test_coal coalescence_golovin coalescence_hall_davis_no_waals lgrngn_cond lgrngn_cond_substepping lgrngn_cond_substepping_test lgrngn_cond_substepping_plot puddle) # the coalescence_onishi_hall test is turned off as not specific enough
  #TODO: indicate that tests depend on the lib
  add_test(
    NAME ${test}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/python" 
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/python/physics/${test}.py"
  )
endforeach()


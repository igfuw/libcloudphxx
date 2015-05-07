set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

#TODO: add test_cuda
add_subdirectory(theta)
add_subdirectory(particles)
add_subdirectory(fake_units)
add_subdirectory(bisect)
add_subdirectory(zip)
add_subdirectory(python)
add_subdirectory(common)
add_subdirectory(blk2m_hello_world)

# TODO: target_compile_options() // added to CMake on Jun 3rd 2013

set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

add_subdirectory(parcel)
add_subdirectory(theta)
add_subdirectory(particles)
add_subdirectory(fake_units)
add_subdirectory(bisect)
add_subdirectory(zip)

# TODO: target_compile_options() // added to CMake on Jun 3rd 2013

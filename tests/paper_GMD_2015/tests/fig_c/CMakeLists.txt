find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_BINARY_DIR})
target_link_libraries(plot ${Boost_LIBRARIES})
target_link_libraries(plot ${HDF5_LIBRARIES})
target_include_directories(plot PUBLIC ${HDF5_INCLUDE_DIRS})


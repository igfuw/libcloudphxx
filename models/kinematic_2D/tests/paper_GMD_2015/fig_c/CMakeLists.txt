add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_BINARY_DIR})
target_link_libraries(plot ${libmpdataxx_LIBRARIES})
target_include_directories(plot PUBLIC ${libmpdataxx_INCLUDE_DIRS})


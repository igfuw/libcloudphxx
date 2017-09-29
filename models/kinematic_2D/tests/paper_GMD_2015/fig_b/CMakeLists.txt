add_executable(calc_b calc.cpp)
target_link_libraries(calc_b ${libmpdataxx_LIBRARIES})
target_include_directories(calc_b PUBLIC ${libmpdataxx_INCLUDE_DIRS})
add_test(calc_b calc_b ${CMAKE_BINARY_DIR})

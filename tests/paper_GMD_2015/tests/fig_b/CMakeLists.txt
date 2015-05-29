add_executable(calc_b calc.cpp)
add_test(calc_b calc_b ${CMAKE_BINARY_DIR})

find_package(Boost COMPONENTS system timer filesystem REQUIRED)
target_link_libraries(calc_b ${Boost_LIBRARIES})

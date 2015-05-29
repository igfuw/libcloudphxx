find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_BINARY_DIR})
target_link_libraries(plot ${Boost_LIBRARIES})

add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_CURRENT_BINARY_DIR})



find_package(Boost COMPONENTS iostreams system REQUIRED)
target_link_libraries(plot ${Boost_LIBRARIES})
target_include_directories(plot PUBLIC ${Boost_INCLUDE_DIRS})
#else()
#    message(FATAL_ERROR "Boost (or some of its components) not found.")
#endif()
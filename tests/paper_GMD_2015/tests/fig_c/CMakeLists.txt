find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_BINARY_DIR})
target_link_libraries(plot ${Boost_LIBRARIES})

# if using local version of libcloudphxx use local headers
if (libcloudphxx_local)
  target_include_directories(plot PUBLIC ${CMAKE_SOURCE_DIR}/../../include)
endif()

find_package(Boost COMPONENTS iostreams system REQUIRED)

add_executable(plot plot.cpp)
add_test(plot plot ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(plot ${Boost_LIBRARIES})
target_include_directories(plot PUBLIC ${Boost_INCLUDE_DIRS})


add_executable(parcel parcel.cpp)
target_link_libraries(parcel PRIVATE clphxx::cloudphxx_lgrngn)
target_compile_features(parcel PUBLIC cxx_std_20)
target_link_libraries(parcel PRIVATE ${Boost_LIBRARIES})
target_include_directories(parcel PUBLIC ${Boost_INCLUDE_DIRS})

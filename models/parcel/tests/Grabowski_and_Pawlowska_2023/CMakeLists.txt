find_package(Boost COMPONENTS program_options REQUIRED)

add_executable(parcel parcel.cpp)
target_link_libraries(parcel PRIVATE clphxx::cloudphxx_lgrngn)
target_compile_features(parcel PUBLIC cxx_std_20)
target_link_libraries(parcel PRIVATE ${Boost_LIBRARIES})
target_include_directories(parcel PUBLIC ${Boost_INCLUDE_DIRS})

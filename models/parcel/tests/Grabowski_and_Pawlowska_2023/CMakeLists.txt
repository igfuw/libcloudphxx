add_executable(parcel parcel.cpp)
target_compile_features(parcel PUBLIC cxx_std_20)

target_link_libraries(parcel PRIVATE clphxx::cloudphxx_lgrngn)

find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(parcel PRIVATE ${Boost_LIBRARIES})
target_include_directories(parcel PUBLIC ${Boost_INCLUDE_DIRS})

find_library(netcdfcxx_LIBRARIES netcdf_c++4 REQUIRED)
target_link_libraries(parcel PRIVATE ${netcdfcxx_LIBRARIES})
find_path(netcdfxx_INCLUDE_DIRS netcdf)
target_include_directories(parcel PUBLIC ${netcdfxx_INCLUDE_DIRS})
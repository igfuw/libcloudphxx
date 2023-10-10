add_executable(parcel parcel.cpp)
target_compile_features(parcel PUBLIC cxx_std_20)

target_link_libraries(parcel PRIVATE clphxx::cloudphxx_lgrngn)

find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(parcel PRIVATE ${Boost_LIBRARIES})
target_include_directories(parcel PUBLIC ${Boost_INCLUDE_DIRS})

#NOTE: on Debian, netcdf_c++4 may not link netcdf libs (https://github.com/Unidata/netcdf-cxx4/issues/86)
#      to fix parcel linkage, add -lnetcdf to CMakeFiles/parcel.dir/link.txt in Grabowski_and_Pawlowska_2023 build directory
find_library(netcdfcxx_LIBRARIES netcdf_c++4 REQUIRED)
target_link_libraries(parcel PRIVATE ${netcdfcxx_LIBRARIES})
find_path(netcdfxx_INCLUDE_DIRS netcdf)
target_include_directories(parcel PUBLIC ${netcdfxx_INCLUDE_DIRS})

find_package(Thrust REQUIRED)
thrust_create_target(Thrust)
target_link_libraries(parcel PRIVATE Thrust)
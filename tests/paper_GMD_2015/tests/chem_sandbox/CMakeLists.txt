add_executable(calc_chem calc_chem.cpp)
add_test(calc_chem calc_chem ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(micro lgrngn lgrngn_spec lgrngn_chem)
  add_executable(plot_${micro} plot_${micro}.cpp)
  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(plot_${micro} ${Boost_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${Boost_INCLUDE_DIRS})
endforeach()

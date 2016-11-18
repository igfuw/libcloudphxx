add_executable(calc_chem calc_chem.cpp)
add_test(calc_chem calc_chem ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(case lgrngn_ch lgrngn_ch_spec lgrngn_ch_chem)
  add_executable(plot_${case} plot_${case}.cpp)
  add_test(plot_${case} plot_${case} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${case} ${HDF5_LIBRARIES})
  target_include_directories(plot_${case} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(plot_${case} ${Boost_LIBRARIES})
  target_include_directories(plot_${case} PUBLIC ${Boost_INCLUDE_DIRS})
endforeach()

add_executable(calc_chem calc_chem.cpp)

foreach(case lgrngn_ch lgrngn_ch_chem)
  add_executable(plot_${case} plot_${case}.cpp)
  target_link_libraries(plot_${case} ${libmpdataxx_LIBRARIES})
  target_include_directories(plot_${case} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
endforeach()

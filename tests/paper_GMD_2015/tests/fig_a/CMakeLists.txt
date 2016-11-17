# tests from GMD_2015 paper (should be run on CUDA)
add_executable(calc calc.cpp)
add_test(calc calc ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(micro blk_1m blk_2m lgrngn)
  add_executable(plot_${micro} plot_${micro}.cpp)
  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(plot_${micro} ${Boost_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${Boost_INCLUDE_DIRS})
endforeach()

add_executable(plot_lgrngn_spec plot_lgrngn_spec.cpp)
add_test(plot_lgrngn_spec plot_lgrngn_spec ${CMAKE_BINARY_DIR})
target_link_libraries(plot_lgrngn_spec ${HDF5_LIBRARIES})
target_include_directories(plot_lgrngn_spec PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(plot_lgrngn_spec ${Boost_LIBRARIES})
target_include_directories(plot_lgrngn_spec PUBLIC ${Boost_INCLUDE_DIRS})

# tests for travis (or any other case where libcloud can't use CUDA)
add_executable(travis_calc_blk travis_calc_blk.cpp)
add_test(travis_calc_blk travis_calc_blk ${CMAKE_BINARY_DIR})

add_executable(travis_calc_lgrngn travis_calc_lgrngn.cpp)
add_test(travis_calc_lgrngn travis_calc_lgrngn ${CMAKE_BINARY_DIR})

add_test(travis_2D_kin_cloud_dims_blk bash -c "
  for dir in out_blk_1m out_blk_2m; do
    echo 'comparing const.h5'
    if h5diff $dir/const.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5 | grep 'not comparable'; then
      echo 'the subsequent diff test will not detect differences -- Fail'
      exit 1      
    else          # If h5diff encounters files that cannot be compared (for example with different dimesions)
      echo 'OK'
      exit 0      # it does not return exit 1. Instead it prints out a warning that 'Some files are not comparable'
    fi
  done
")

add_test(travis_2D_kin_cloud_diff_blk_1m bash -c "
  for dir in out_blk_1m; do
    echo 'comparing const.h5'                                                                                                &&
    h5diff --delta=2e-7  -v2 $dir/const.h5           ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                       &&
    echo 'comparing timestep0000000000.h5'                                                                                   &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rv  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rr  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rc  &&
    h5diff --delta=0.1   -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /th  &&
    echo 'comparing timestep0000009000.h5'                                                                                   &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /rv  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /rr  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /rc  &&
    h5diff --delta=0.1   -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /th  || exit 1;
  done 
")

add_test(travis_2D_kin_cloud_diff_blk_2m bash -c "
  for dir in out_blk_2m; do
    echo 'comparing const.h5'                                                                                                   &&
    h5diff --delta=2e-7     -v2 $dir/const.h5              ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                    &&
    echo 'comparing timestep0000000000.h5'                                                                                      &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rv  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rr  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /rc  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /nr  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /nc  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  /th  &&
    echo 'comparing timestep0000009000.h5'                                                                                      &&
    h5diff --relative=0.02  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /rv  &&
    h5diff --relative=0.02  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /rr  &&
    h5diff --relative=0.02  -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /nr  &&
    h5diff --delta=0.4      -v2 $dir/timestep0000009000.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000009000.h5  /th  || exit 1;
  done  
")

add_test(travis_2D_kin_cloud_dims_lgrngn bash -c "
    echo 'comparing dimensions in const.h5'
    if h5diff out_lgrngn/const.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/out_lgrngn/travis_const.h5 | grep 'not comparable'; then
      echo 'the subsequent diff test will not detect differences -- Fail'
      exit 1      
    else          # If h5diff encounters files that cannot be compared (for example with different dimesions)
      echo 'OK'
      exit 0      # it does not return exit 1. Instead it prints out a warning that 'Some files are not comparable'.
    fi
")

add_test(travis_2D_kin_cloud_diff_lgrngn bash -c "
  for dir in out_lgrngn; do
    echo 'comparing const.h5'                                                                                                                && 
    h5diff --delta=2e-7     -v2 $dir/const.h5                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_const.h5                       &&
    echo 'comparing timestep0000000000.h5'                                                                                                   &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5    ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000000.h5 /th      &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5    ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000000.h5 /rv      &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5    ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000000.h5 /sd_conc &&
    echo 'comparing timestep0000000020.h5'                                                                                                   &&
    h5diff --relative=0.001 -v2 $dir/timestep0000000020.h5   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000020.h5 /th       &&
    h5diff --relative=0.001 -v2 $dir/timestep0000000020.h5   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000020.h5 /rv       &&
    h5diff --delta=1        -v2 $dir/timestep0000000020.h5   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/travis_timestep0000000020.h5 /sd_conc || exit 1;
  done  
")

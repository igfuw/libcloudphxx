# CMake generated Testfile for 
# Source directory: /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a
# Build directory: /home/piotr/praca/libcloudphxx/models/kinematic_2D/build/tests/paper_GMD_2015/fig_a
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(calc "calc" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(plot_blk_1m "plot_blk_1m" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(plot_blk_2m "plot_blk_2m" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(plot_lgrngn "plot_lgrngn" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(plot_lgrngn_spec "plot_lgrngn_spec" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(travis_calc_blk "travis_calc_blk" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(travis_calc_lgrngn "travis_calc_lgrngn" "/home/piotr/praca/libcloudphxx/models/kinematic_2D/build")
add_test(travis_2D_kin_cloud_dims_blk "bash" "-c" "
  for dir in out_blk_1m out_blk_2m; do
    echo 'comparing const.h5'
    if h5diff $dir/const.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/const.h5 | grep 'not comparable'; then
      echo 'the subsequent diff test will not detect differences -- Fail'
      exit 1      
    else          # If h5diff encounters files that cannot be compared (for example with different dimesions)
      echo 'OK'
      exit 0      # it does not return exit 1. Instead it prints out a warning that 'Some files are not comparable'
    fi
  done
")
add_test(travis_2D_kin_cloud_diff_blk_1m "bash" "-c" "
  for dir in out_blk_1m; do
    echo 'comparing const.h5'                                                                                                &&
    h5diff --delta=2e-7  -v2 $dir/const.h5           /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/const.h5                       &&
    echo 'comparing timestep0000000000.h5'                                                                                   &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rv  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rr  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rc  &&
    h5diff --delta=0.1   -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /th  &&
    echo 'comparing timestep0000009000.h5'                                                                                   &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rv  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rr  &&
    h5diff --delta=2e-5  -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rc  &&
    h5diff --delta=0.1   -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /th  || exit 1;
  done 
")
add_test(travis_2D_kin_cloud_diff_blk_2m "bash" "-c" "
  for dir in out_blk_2m; do
    echo 'comparing const.h5'
    h5diff --delta=2e-7     -v2 $dir/const.h5              /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/const.h5                    &&
    echo 'comparing timestep0000000000.h5'                                                                                     &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rv  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rr  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /rc  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /nr  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /nc  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000000000.h5  /th  &&
    echo 'comparing timestep0000009000.h5'                                                                                     &&
    h5diff --relative=0.02  -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rv  &&
    h5diff --delta=12e-6     -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rr  &&
    h5diff --delta=4e-6     -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /rc  &&
    h5diff --delta=0.4      -v2 $dir/timestep0000009000.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/timestep0000009000.h5  /th  || exit 1;
  done  
")
add_test(travis_2D_kin_cloud_dims_lgrngn "bash" "-c" "
    echo 'comparing dimensions in const.h5'
    if h5diff out_lgrngn/const.h5 /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/out_lgrngn/travis_const.h5 | grep 'not comparable'; then
      echo 'the subsequent diff test will not detect differences -- Fail'
      exit 1      
    else          # If h5diff encounters files that cannot be compared (for example with different dimesions)
      echo 'OK'
      exit 0      # it does not return exit 1. Instead it prints out a warning that 'Some files are not comparable'.
    fi
")
add_test(travis_2D_kin_cloud_diff_lgrngn "bash" "-c" "
  for dir in out_lgrngn; do
    echo 'comparing const.h5'
    h5diff --delta=1e-6     -v2 $dir/const.h5                /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_const.h5                       &&
    echo 'comparing timestep0000000000.h5'                                                                                                  &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000000.h5 /th      &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000000.h5 /rv      &&
    h5diff --relative=1e-9  -v2 $dir/timestep0000000000.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000000.h5 /sd_conc &&
    echo 'comparing timestep0000000020.h5'                                                                                                  &&
    h5diff --relative=0.001 -v2 $dir/timestep0000000020.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000020.h5 /th      &&
    h5diff --relative=0.002 -v2 $dir/timestep0000000020.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000020.h5 /rv      &&
    h5diff --delta=1        -v2 $dir/timestep0000000020.h5   /home/piotr/praca/libcloudphxx/models/kinematic_2D/tests/paper_GMD_2015/fig_a/refdata/$dir/travis_timestep0000000020.h5 /sd_conc || exit 1;
  done  
")

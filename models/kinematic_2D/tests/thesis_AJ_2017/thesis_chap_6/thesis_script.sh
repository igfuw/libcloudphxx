#!/bin/bash

# Runs all the simulations and plots the results from chapter 6 of PhD thesis of Jaruga.
# The simulations take up to 12 hours (on cuda-k-2 server using CUDA and one thread)
# It might be useful to run them in the background: nohup ./thesis_script.sh &

cd ../../../build/tests/thesis_AJ_2017/thesis_chap_6/

export OMP_NUM_THREADS=1

## run all the thesis test cases
./calc_chem ../../../

## plot the quicklook plots
./plot_lgrngn_ch ../../
./plot_lgrngn_ch_chem ../../

# plot the python scripts for pH and size distributions
python ../../../../tests/thesis_AJ_2017/thesis_chap_6/chem_plots/ph_plot.py
python ../../../../tests/thesis_AJ_2017/thesis_chap_6/chem_plots/rain_histograms_single.py
python ../../../../tests/thesis_AJ_2017/thesis_chap_6/chem_plots/rain_histograms_all.py

# check how much H2SO4 mass was created and what is the error
python ../../../../tests/thesis_AJ_2017/thesis_chap_6/chem_plots/how_effective.py
python ../../../../tests/thesis_AJ_2017/thesis_chap_6/chem_plots/test_mole_const.py

# thesis_plots folder will contain all plots included in thesis
mkdir -p thesis_plots

# copy quicklook plots from case1 and case4
for sim_run in case1 case4; do
  for el in nc na ef rd rr sd_conc; do
    dir_name=out_${sim_run}.plot/${el}/;
    mkdir -p thesis_plots/${dir_name};
    for sim_time in 10000 11800; do
      cp ${dir_name}/${sim_time}.svg thesis_plots/${dir_name}/${sim_time}.svg;
    done
  done
done

# copy quicklook plots from the base case
for el in nc na ef rd rr sd_conc SO2g O3g H2O2g S_VI_aq H2O2_aq O3_aq; do
  dir_name=out_case_base.plot/${el}/;
  mkdir -p thesis_plots/${dir_name};
  for sim_time in 10000 11800; do
    cp ${dir_name}/${sim_time}.svg thesis_plots/${dir_name}/${sim_time}.svg;
  done
done

# copy quicklook plots from case2
for kernel in hall hall_davis_no_waals hall_pinsky_stratocumulus onishi_hall onishi_hall_davis_no_waals vohl_davis_no_waals; do
  dir_name=out_case2_${kernel}_42.plot/rr/;
  mkdir -p thesis_plots/${dir_name};
  cp ${dir_name}/11800.svg thesis_plots/${dir_name}/11800.svg;
done

# copy pH and size distribution plots
cp -r plots_of_pH/ thesis_plots/.
cp -r plots_of_size_distr/ thesis_plots/.

# pack all plots included in thesis
tar -cvf thesis_plots.tar thesis_plots

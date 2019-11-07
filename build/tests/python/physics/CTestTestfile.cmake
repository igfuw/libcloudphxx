# CMake generated Testfile for 
# Source directory: /mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics
# Build directory: /mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/tests/python/physics
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_coal "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/test_coal.py")
set_tests_properties(test_coal PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")
add_test(coalescence_golovin "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/coalescence_golovin.py")
set_tests_properties(coalescence_golovin PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")
add_test(coalescence_hall_davis_no_waals "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/coalescence_hall_davis_no_waals.py")
set_tests_properties(coalescence_hall_davis_no_waals PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")
add_test(lgrngn_cond "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/lgrngn_cond.py")
set_tests_properties(lgrngn_cond PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")
add_test(lgrngn_cond_substepping "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/lgrngn_cond_substepping.py")
set_tests_properties(lgrngn_cond_substepping PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")
add_test(puddle "/usr/bin/python2.7" "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/tests/python/physics/puddle.py")
set_tests_properties(puddle PROPERTIES  WORKING_DIRECTORY "/mnt/local/pzmij/UWLCM_singu/singularity/libcloudphxx/build/bindings/python")

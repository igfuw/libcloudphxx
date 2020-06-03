add_executable(mpi_adve_test mpi_adve_test.cpp)
target_link_libraries(mpi_adve_test cloudphxx_lgrngn)

add_test(NAME mpi_adve_test_np1 COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 mpirun -np 1  ./mpi_adve_test  -c 0 -d 1)
add_test(NAME mpi_adve_test_np2 COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 mpirun -np 2  ./mpi_adve_test  -c 0 -d 1)
add_test(NAME mpi_adve_test_np3 COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 mpirun -np 3  ./mpi_adve_test  -c 0 -d 1)
add_test(NAME mpi_adve_test_np4 COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 mpirun -np 4  ./mpi_adve_test  -c 0 -d 1)

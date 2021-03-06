function(mpi_test TESTNAME PROCS EXE)
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${EXE} ${ARGN}
  )
endfunction(mpi_test)

mpi_test(adjacencies_2D 1  ./adjacencies ${CMAKE_SOURCE_DIR}/meshes/tri8.osh)
mpi_test(tags_2D 1  ./tags ${CMAKE_SOURCE_DIR}/meshes/tri8.osh)
mpi_test(interpolation_3D 1  ./interpolation ${CMAKE_SOURCE_DIR}/meshes/cube7k.osh)
mpi_test(classification_3D 1  ./classification ${CMAKE_SOURCE_DIR}/meshes/cube7k.osh)
mpi_test(reduction_2D 4  ./reduction ${CMAKE_SOURCE_DIR}/meshes/cube7k_4p.osh)
mpi_test(synchronization_2D 4  ./synchronization  ${CMAKE_SOURCE_DIR}/meshes/tri8_4p.osh)
mpi_test(partitioning_2D 4  ./partitioning  ${CMAKE_SOURCE_DIR}/meshes/tri8_4p.osh)
mpi_test(ghosting_2D 4  ./ghosting  ${CMAKE_SOURCE_DIR}/meshes/cube7k_4p.osh)
mpi_test(hypercube_3D 1  ./hypercube  ${CMAKE_SOURCE_DIR}/meshes/hypercube.osh)

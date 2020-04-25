#pragma once

namespace
{
  bool ran_with_mpi()
  {
   return (
     // mpich
     std::getenv("PMI_RANK") != NULL ||
     // openmpi
     std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
     // lam
     std::getenv("LAMRANK") != NULL ||
     // mvapich2
     std::getenv("MV2_COMM_WORLD_RANK") != NULL
   );
  }
};

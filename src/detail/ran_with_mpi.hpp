#pragma once

namespace
{
  bool ran_with_mpi()
  {
   return (
     // mpich
     std::getenv("PMI_RANK") != NULL && std::atoi(std::getenv("PMI_SIZE")) > 1 ||
     // openmpi
     std::getenv("OMPI_COMM_WORLD_RANK") != NULL && std::atoi(std::getenv("OMPI_COMM_WORLD_SIZE")) > 1 ||
     // lam
     std::getenv("LAMRANK") != NULL && std::atoi(std::getenv("LAMSIZE")) > 1 ||
     // mvapich2
     std::getenv("MV2_COMM_WORLD_RANK") != NULL && std::atoi(std::getenv("MV2_COMM_WORLD_SIZE")) > 1
   );
  }
};

#pragma once

// Error handling macro
#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    std::cerr << "MPI error calling \""#call"\"\n"; \
    MPI_Abort(MPI_COMM_LIBCLOUD, -1);}

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      namespace
      {
        MPI_Comm MPI_COMM_LIBCLOUD; // our copy of MPI_COMM_WORLD to have our separate context
        bool MPI_COMM_LIBCLOUD_INITIALIZED = false;

        int mpi_init(const int &req_tlvl)
        {
          int initialized, prov_tlvl;

          MPI_CHECK(MPI_Initialized(&initialized));
          if(!initialized)
          {
            MPI_CHECK(MPI_Init_thread(NULL, NULL, req_tlvl, &prov_tlvl));
          }
          else
          {
            MPI_CHECK(MPI_Query_thread(&prov_tlvl));
          }

          if(!MPI_COMM_LIBCLOUD_INITIALIZED)
          {
            MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_LIBCLOUD);
            MPI_COMM_LIBCLOUD_INITIALIZED = true;
          }

          return prov_tlvl;
        }

        int mpi_init(const int &req_tlvl, int &rank, int &size)
        {
          int prov_tl = mpi_init(req_tlvl);
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_LIBCLOUD, &rank));
          MPI_CHECK(MPI_Comm_size(MPI_COMM_LIBCLOUD, &size));
          return prov_tl;
        }
      }
    };
  };
};


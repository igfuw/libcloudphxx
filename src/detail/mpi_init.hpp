#pragma once

// Error handling macro
#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    std::cerr << "MPI error calling \""#call"\"\n"; \
    MPI_Abort(MPI_COMM_WORLD, -1);}

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      namespace
      {
        int mpi_init(const int &req_tlvl)
        {
          // TODO: create a communicator for 'our' processes?
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
          return prov_tlvl;
        }

        int mpi_init(const int &req_tlvl, int &rank, int &size)
        {
          int prov_tl = mpi_init(req_tlvl);
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
          MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
          return prov_tl;
        }
      }
    };
  };
};


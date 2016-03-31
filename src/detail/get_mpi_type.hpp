#pragma once 

#if defined(USE_MPI)
namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      namespace
      {
        // mpi message tags
        enum {tag_n_lft = 654321, tag_real_lft, tag_n_rgt, tag_real_rgt}; // hope other libs dont use these tags, TODO: using separate communicator would help?

        template<typename real_t>
        MPI_Datatype get_mpi_type()
        {
          throw std::runtime_error("Unsupported MPI datatype");
        }

        template<>
        MPI_Datatype get_mpi_type<float>()
        {
          return MPI_FLOAT;
        }

        template<>
        MPI_Datatype get_mpi_type<double>()
        {
          return MPI_DOUBLE;
        }

        template<>
        MPI_Datatype get_mpi_type<unsigned long long>()
        {
          return MPI_UNSIGNED_LONG_LONG;
        }
      }
    }
  }
}
#endif


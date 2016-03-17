#pragma once

// Error handling macros
#if defined(USE_MPI)
  #define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
      std::cerr << "MPI error calling \""#call"\"\n"; \
      MPI_Abort(MPI_COMM_WORLD, -1);}
#endif



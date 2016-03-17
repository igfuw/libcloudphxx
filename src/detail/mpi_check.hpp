#pragma once

// Error handling macros
#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    std::cerr << "MPI error calling \""#call"\"\n"; \
    MPI_Abort(MPI_COMM_WORLD, -1);}



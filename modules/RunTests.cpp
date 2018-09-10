#include <iostream>

// include the google testing header
#include <gtest/gtest.h>

#include "MUQ/config.h"

#if MUQ_HAS_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

#if MUQ_HAS_MPI
  const int ierr = MPI_Init(nullptr, nullptr);
#endif

  const int res = RUN_ALL_TESTS();   
  
#if MUQ_HAS_MPI
  MPI_Finalize();
#endif
  
  return res;
}

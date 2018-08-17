#include <iostream>

// include the google testing header
#include <gtest/gtest.h>

#include "MUQ/config.h"

#if MUQ_HAS_MPI==0
#error
#endif

#include <mpi.h>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  const int ierr = MPI_Init(nullptr, nullptr);

  const int res = RUN_ALL_TESTS();   
  
  MPI_Finalize();
  
  return res;
}

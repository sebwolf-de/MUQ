#include <iostream>

// include the google testing header
#include <gtest/gtest.h>

#include "MUQ/config.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

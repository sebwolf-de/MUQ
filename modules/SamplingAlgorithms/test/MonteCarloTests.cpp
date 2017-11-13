#include <gtest/gtest.h>

#include "MUQ/SamplingAlgorithms/MonteCarlo.h"

using namespace muq::SamplingAlgorithms;

TEST(MonteCarlo, Setup) {
  // create an instance of Monte Carlo
  auto mc = std::make_shared<MonteCarlo>();
}

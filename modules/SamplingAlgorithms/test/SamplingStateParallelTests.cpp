#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>

#include <cereal/access.hpp>

#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::SamplingAlgorithms;

TEST(SamplingState, Serialization) {
  // random state
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(3);
  //auto state = std::make_shared<SamplingState>(vec, 1.0);
  SamplingState state(vec, 1.0);

  {
    //cereal::BinaryOutputArchive oarchive(std::cout);
    //oarchive(state);

    //state->serialize<cereal::XMLInputArchive>(iarchive);
  }

}

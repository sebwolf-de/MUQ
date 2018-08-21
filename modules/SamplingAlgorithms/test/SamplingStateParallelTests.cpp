#include "MUQ/config.h"

#if MUQ_HAS_PARCER

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <fstream>
#include <iostream>

#include <cereal/archives/binary.hpp>

#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::SamplingAlgorithms;

TEST(SamplingState, Serialization) {

  const Eigen::VectorXd vec = Eigen::VectorXd::Random(3);
  SamplingState state(vec, 0.5);

  state.meta["Something"] = double(1.5);
  state.meta["SomeString"] = std::string("BlahBlahBlah");

  {
    std::ofstream out("eigen.cereal", std::ios::binary);
    cereal::BinaryOutputArchive archive_o(out);
    archive_o(state);
  }

  SamplingState state2;
  {
    std::ifstream in("eigen.cereal", std::ios::binary);
    cereal::BinaryInputArchive archive_i(in);
    archive_i(state2);
  }

  for(int i=0; i<vec.size(); ++i)
    EXPECT_DOUBLE_EQ(state.state.at(0)(i), state2.state.at(0)(i));

  EXPECT_DOUBLE_EQ(state.weight, state2.weight);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(state.meta["Something"]), boost::any_cast<double>(state2.meta["Something"]));
  EXPECT_EQ(boost::any_cast<std::string>(state.meta["SomeString"]), boost::any_cast<std::string>(state2.meta["SomeString"]));
}


#endif

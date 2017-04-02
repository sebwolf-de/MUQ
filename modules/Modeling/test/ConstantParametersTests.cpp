#include "gtest/gtest.h"

#include "MUQ/Modeling/ConstantParameters.h"

using namespace muq::Modeling;

/// A generic object for testing purposes
struct AnObject {
  AnObject(double const value) : value(value) {}

  const double value;
};

TEST(ConstantParameters, BasicTest) {
  // create an object
  auto obj = std::make_shared<AnObject>(2.0);

  // the outputs to the constant parameter
  std::vector<boost::any> outputs({1.0, (std::string)"string", obj});

  // create a constant parameter
  auto para = std::make_shared<ConstantParameters>(outputs);

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, 3);

  // evaluate the outputs
  auto outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 3);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs[0]), 1.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs[1]).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outs[2])->value, 2.0);
}

TEST(ConstantParameters, RecursiveConstructor) {
  // create an object
  auto obj = std::make_shared<AnObject>(2.0);

  // create a constant parameter
  auto para = std::make_shared<ConstantParameters>(1.0, (std::string)"string", obj);

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, 3);

  // evaluate the outputs
  auto outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 3);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs[0]), 1.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs[1]).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outs[2])->value, 2.0);
}

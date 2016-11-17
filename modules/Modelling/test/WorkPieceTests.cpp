#include "gtest/gtest.h"

#include "MUQ/Modelling/WorkPiece.h"

//WorkParent(WrkPrnt, double)

WorkParent(TestIntParent, 3, int, double, std::string, double, int)

class TestIntMod : public muq::Modelling::TestIntParent {
public:
  /**
     @param[in] numIns The number of inputs
   */
  TestIntMod() : TestIntParent() {}

  virtual ~TestIntMod() {}
private:
};

TEST(WorkPiece, Build) {
  // create the test WorkPiece
  auto test = std::make_shared<TestIntMod>();

  test->Evaluate("a string", 2.0, 3);
  //std::cout << std::get<0>(test->Evaluate()) << std::endl;
  
  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->NumInputs(), 3);

  // make sure the number of outputs matches what we expect
  EXPECT_EQ(test->NumOutputs(), 2);

  // make sure the input types are correct
  const std::vector<std::string> ins = test->InputTypes();
  EXPECT_TRUE(ins.at(0).compare("std::string")==0);
  EXPECT_TRUE(ins.at(1).compare("double")==0);
  EXPECT_TRUE(ins.at(2).compare("int")==0);

  // make sure the output types are correct
  const std::vector<std::string> outs = test->OutputTypes();
  EXPECT_TRUE(outs.at(0).compare("int")==0);
  EXPECT_TRUE(outs.at(1).compare("double")==0);
}

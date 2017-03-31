#include "gtest/gtest.h"

#include "WorkPieceTestClasses.h"

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class WorkPieceTests : public::testing::Test {
public:

  /// Default constructor
  WorkPieceTests() {}

  /// Default destructor
  virtual ~WorkPieceTests() {}

  /// A string type object to input
  const std::string s = "a string";

  /// A double type object to input
  const double a = 2.0;

  /// An unsigned int type object to input
  const unsigned int b = 3;

private:
};

TEST_F(WorkPieceTests, UnfixedInputOutputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<UnfixedMod>();

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // evaluate with zero inputs and three outputs
  std::vector<boost::any> outputs = test->Evaluate();

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 3);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare("hello!")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 3.0);
  EXPECT_EQ(boost::any_cast<int>(outputs[2]), 6);

  // evaluate with one input and one output
  outputs = test->Evaluate(s);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  
  // evaluate with two inputs and zero outputs
  outputs = test->Evaluate(a, b);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 0);
}

TEST_F(WorkPieceTests, FixedInputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedInsMod>(3);

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // evaluate with a bool that will return both outputs
  std::vector<boost::any> outputs = test->Evaluate(s, b, true);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_EQ(boost::any_cast<unsigned int>(outputs[1]), b);

  // evaluate with a bool that will return one output
  outputs = test->Evaluate(s, b, false);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_EQ(boost::any_cast<unsigned int>(outputs[0]), b);
}

TEST_F(WorkPieceTests, FixedOutputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedOutsMod>(2);

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // evaluate with no inputs
  std::vector<boost::any> outputs = test->Evaluate();

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), 1.0);
  EXPECT_EQ(boost::any_cast<int>(outputs[1]), 2);

  // evaluate with one inputs
  outputs = test->Evaluate(3);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), 1.0);
  EXPECT_EQ(boost::any_cast<int>(outputs[1]), 3);
}

TEST_F(WorkPieceTests, FixedInputOutputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(3, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);
}

TEST_F(WorkPieceTests, FixedInputTypes) {
  // the input types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInTypeMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
}

TEST_F(WorkPieceTests, FixedOutputTypes) {
  // the output types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedOutTypeMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);
}

TEST_F(WorkPieceTests, SomeFixedInputTypes) {
  // the input types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[2] = typeid(std::shared_ptr<AnObject>).name();

  // create the test WorkPiece
  auto test = std::make_shared<SomeFixedInTypeMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj, 9.0, 'c');

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
}

TEST_F(WorkPieceTests, SomeFixedOutputTypes) {
  // the output types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<SomeFixedOutTypeMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a, 'c');

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 3);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);
  EXPECT_TRUE(boost::any_cast<bool>(outputs[2]));

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s);
  
  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
}

TEST_F(WorkPieceTests, FixedInputTypesOutputNum) {
  // the input types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInTypeOutNumMod>(types, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  const std::string outString = "second string";
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[1]).compare(outString)==0);
}

TEST_F(WorkPieceTests, FixedOutputTypesInputNum) {
  // the output types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedOutTypeInNumMod>(types, 3);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, obj, 'a'); // use a random character as the third input (fixed number, not fixed type)

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);
}

TEST_F(WorkPieceTests, FixedTypes) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(std::shared_ptr<AnObject>).name(), typeid(double).name()});

  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedTypesMod>(inTypes, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, obj, 1.0);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);
}

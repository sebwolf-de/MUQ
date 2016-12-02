#include "gtest/gtest.h"

#include "MUQ/Modelling/WorkPiece.h"

struct AnObject {
  AnObject(double const b) : value(b) {}

  bool flag; 

  const double value;
};

class UnfixedMod : public muq::Modelling::WorkPiece {
public:

  UnfixedMod() : WorkPiece() {}

  virtual ~UnfixedMod() {}

private:

  virtual void EvaluateImpl() override {
    switch( inputs.size() ) {
    case 0 : {
      const std::string hi = "hello!";

      outputs.resize(3);

      outputs[0] = hi;
      outputs[1] = 3.0;
      outputs[2] = 6;

      return;
    } case 1: {
      outputs.resize(1);
      outputs[0] = boost::any_cast<std::string>(inputs[0]);

      return;
    } default:
      return;
    }
  }    
};

class FixedInsMod : public muq::Modelling::WorkPiece {
public:

  FixedInsMod(unsigned int numIns) : WorkPiece(numIns) {}

  virtual ~FixedInsMod() {}

private:

  virtual void EvaluateImpl() override {
    if( boost::any_cast<bool>(inputs[2]) ) {
      outputs.resize(2);

      outputs[0] = boost::any_cast<std::string>(inputs[0]);
      outputs[1] = boost::any_cast<unsigned int>(inputs[1]);
    } else {
      outputs.resize(1);

      outputs[0] = boost::any_cast<unsigned int>(inputs[1]);
    }
  }    
};

class FixedOutsMod : public muq::Modelling::WorkPiece {
public:

  FixedOutsMod(unsigned int numIns) : WorkPiece(numIns, false) {}

  virtual ~FixedOutsMod() {}

private:

  virtual void EvaluateImpl() override {
    outputs.resize(numOutputs);
    
    outputs[0] = 1.0;
    outputs[1] = inputs.size()>0? boost::any_cast<int>(inputs[0]) : 2;
  }    
};

class FixedInOutMod : public muq::Modelling::WorkPiece {
public:

  FixedInOutMod(unsigned int const numIns, unsigned int numOuts) : WorkPiece(numIns, numOuts) {}

  virtual ~FixedInOutMod() {}

private:

  virtual void EvaluateImpl() override {
    outputs = std::vector<boost::any>();
  }
};

class FixedInTypeMod : public muq::Modelling::WorkPiece {
public:

  FixedInTypeMod(std::vector<std::string> const& types) : WorkPiece(types) {}

  virtual ~FixedInTypeMod() {}

private:

  virtual void EvaluateImpl() override {
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);

    if( obj->flag ) {
      outputs.resize(2);
      outputs[0] = s;
      outputs[1] = obj->value*boost::any_cast<double>(inputs[1]);

      return;
    }

    outputs.resize(1);
    outputs[0] = s;
  }
};

class FixedOutTypeMod : public muq::Modelling::WorkPiece {
public:

  FixedOutTypeMod(std::vector<std::string> const& types) : WorkPiece(types, false) {}

  virtual ~FixedOutTypeMod() {}

private:

  virtual void EvaluateImpl() override {
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[1]);

    outputs.resize(2);
    outputs[0] = s;
    
    if( obj->flag ) {
      outputs[1] = obj->value*boost::any_cast<double>(inputs[2]);
    } else {
      outputs[1] = (double)obj->value;
    }
  }
};

class WorkPieceTests : public::testing::Test {
public:

  WorkPieceTests() {}

  virtual ~WorkPieceTests() {}

  // the inputs
  const std::string s = "a string";
  const double a = 2.0;
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

TEST_F(WorkPieceTests, FixedInputTypesNum) {
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

TEST_F(WorkPieceTests, FixedOutputTypesNum) {
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

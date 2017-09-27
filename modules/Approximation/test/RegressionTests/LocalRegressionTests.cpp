#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

/// A model to approximate (3 input dimensions, 2 output dimensions)
class func : public WorkPiece {
public:

  inline func() :
    WorkPiece(std::vector<std::string>({typeid(Eigen::Vector3d).name()}), // inputs: Eigen::Vector3d
	      std::vector<std::string>({typeid(Eigen::Vector2d).name()})) // outputs: Eigen::Vector2d
  {}

  inline virtual ~func() {}

private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    const Eigen::Vector3d& in = boost::any_cast<Eigen::Vector3d const&>(inputs[0]);
    
    outputs.resize(1);
    outputs[0] = Eigen::Vector2d(in(0)*in(1), in(2));
  }
};

class LocalRegressionTest : public::testing::Test {
public:

  inline LocalRegressionTest() {
    // the function too approximate
    fn = std::make_shared<func>();

    // set the regressor options
    pt::ptree pt;
    pt.put<unsigned int>("LocalRegression.NumNeighbors", 11);
    pt.put<unsigned int>("LocalRegression.Order", 2);

    // create a local regressor
    reg = std::make_shared<LocalRegression>(fn, pt);
  }

  inline virtual ~LocalRegressionTest() {}

  /// The function to approximate
  std::shared_ptr<func> fn;

  /// The local regressor
  std::shared_ptr<LocalRegression> reg;

private:
};

TEST_F(LocalRegressionTest, Basic) {
  // generate some random inputs
  std::vector<Eigen::Vector3d> inputs(25);
  for( auto it=inputs.begin(); it!=inputs.end(); ++it ) { *it = Eigen::Vector3d::Random(); }

  // add the random input points to the cache
  reg->Add(inputs);

  // check the size
  EXPECT_EQ(reg->CacheSize(), inputs.size());

  // the input point
  const Eigen::Vector3d input = Eigen::Vector3d::Random();

  // evaluate the local polynomial approximation
  const std::vector<boost::any>& output = reg->Evaluate(input);
  const Eigen::VectorXd& result = boost::any_cast<Eigen::VectorXd const&>(output[0]);

  // evaluate the truth
  const std::vector<boost::any>& output_truth = fn->Evaluate(input);
  const Eigen::Vector2d& truth = boost::any_cast<Eigen::Vector2d const&>(output_truth[0]);
  
  // the regression and the truth are the same---approximating a quadratic with a quardratic
  EXPECT_NEAR((truth-result).norm(), 0.0, 1.0e-14);
}

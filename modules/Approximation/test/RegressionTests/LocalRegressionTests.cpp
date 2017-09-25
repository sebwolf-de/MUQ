#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/LocalRegression.h"

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
  }
};

class LocalRegressionTest : public::testing::Test {
public:

  inline LocalRegressionTest() {
    // the function too approximate
    fn = std::make_shared<func>();

    // create a local regressor
    reg = std::make_shared<LocalRegression>(fn);
  }

  inline virtual ~LocalRegressionTest() {}

  /// The function to approximate
  std::shared_ptr<func> fn;

  /// The local regressor
  std::shared_ptr<LocalRegression> reg;
  
private:
};

TEST_F(LocalRegressionTest, Basic) {
}

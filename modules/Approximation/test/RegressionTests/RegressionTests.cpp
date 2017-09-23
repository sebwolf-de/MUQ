#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/Regression.h"

using namespace muq::Approximation;

class RegressionTest : public::testing::Test {
public:
  inline RegressionTest() {
    const unsigned int Npts = 100;

    // generate the input points
    ins.resize(100, Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN()));
    for( auto it=ins.begin(); it!=ins.end(); ++it ) { *it = Eigen::Vector2d::Random(); }

    // generate the output points
    outs.resize(ins.size(), Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN()));
    for( unsigned int i=0; i<ins.size(); ++i ) {
      outs[i](0) = ins[i](1)*ins[i](1)*ins[i](1)+ins[i](0)*ins[i](1)-ins[i](0);
      outs[i](1) = ins[i](0)*ins[i](1)*ins[i](1)+ins[i](0)*ins[i](0)+1.5;
    }
  }

  inline virtual ~RegressionTest() {}

  inline virtual void TearDown() override {
    // fit the polynomial coefficients
    reg->Fit<Eigen::Vector2d>(ins, outs);
    
    // points to test the evaluate
    const Eigen::Vector2d x = Eigen::Vector2d::Random();
    const Eigen::Vector2d y = Eigen::Vector2d::Random();
    const Eigen::Vector2d z = Eigen::Vector2d::Random();

    // evaluate the polynomial
    const std::vector<boost::any>& output = reg->Evaluate(x, y, z);
    const Eigen::MatrixXd& result = boost::any_cast<Eigen::MatrixXd const&>(output[0]);
    EXPECT_EQ(result.rows(), 2);    
    EXPECT_EQ(result.cols(), 3);

    // compute the true function values---should be exact (we are using 3rd order to estimate a degree 3 polynomial)
    const Eigen::Vector2d x_true(x(1)*x(1)*x(1)+x(0)*x(1)-x(0), x(0)*x(1)*x(1)+x(0)*x(0)+1.5);
    const Eigen::Vector2d y_true(y(1)*y(1)*y(1)+y(0)*y(1)-y(0), y(0)*y(1)*y(1)+y(0)*y(0)+1.5);
    const Eigen::Vector2d z_true(z(1)*z(1)*z(1)+z(0)*z(1)-z(0), z(0)*z(1)*z(1)+z(0)*z(0)+1.5);

    EXPECT_NEAR((x_true-result.col(0)).norm(), 0.0, 1.0e-12);
    EXPECT_NEAR((y_true-result.col(1)).norm(), 0.0, 1.0e-12);
    EXPECT_NEAR((z_true-result.col(2)).norm(), 0.0, 1.0e-12);
  }

  /// A matrix holding the input points.
  std::vector<Eigen::Vector2d> ins;

  /// A matrix holding the output points.
  std::vector<Eigen::Vector2d> outs;

  /// The order of the polynomial regression
  const unsigned int order = 3;

  /// The object that does the regression
  std::shared_ptr<Regression> reg;
};

TEST_F(RegressionTest, HermiteBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(2, order, Regression::PolynomialBasis::HermiteBasis);
}

TEST_F(RegressionTest, LegendreBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(2, order);
}


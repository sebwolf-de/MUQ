#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/Regression.h"

using namespace muq::Approximation;

class RegressionTest : public::testing::Test {
public:
  inline RegressionTest() {
    const unsigned int Npts = 100;

    // generate the input points
    ins.resize(100, Eigen::VectorXd::Constant(2, std::numeric_limits<double>::quiet_NaN()));
    for( auto it=ins.begin(); it!=ins.end(); ++it ) { *it = Eigen::Vector2d::Random(); }

    // generate the output points
    outs.resize(ins.size());
    for( std::vector<Eigen::VectorXd>::size_type i=0; i<ins.size(); ++i ) { outs[i] = f(ins[i]); }
  }

  inline virtual ~RegressionTest() {}

  inline Eigen::VectorXd f(Eigen::VectorXd const& x) const {
    Eigen::VectorXd out = Eigen::VectorXd::Constant(2, std::numeric_limits<double>::quiet_NaN());
    out(0) = x(1)*x(1)*x(1)+x(0)*x(1)-x(0);
    out(1) = x(0)*x(1)*x(1)+x(0)*x(0)+1.5;

    return out;
  }

  inline virtual void TearDown() override {
    // fit the polynomial coefficients
    reg->Fit(ins, outs);
    
    // points to test the evaluate
    const Eigen::VectorXd x = Eigen::Vector2d::Random();
    const Eigen::VectorXd y = Eigen::Vector2d::Random();
    const Eigen::VectorXd z = Eigen::Vector2d::Random();

    // evaluate the polynomial
    const std::vector<boost::any>& output = reg->Evaluate(x, y, z);
    const Eigen::MatrixXd& result = boost::any_cast<Eigen::MatrixXd const&>(output[0]);
    EXPECT_EQ(result.rows(), 2);    
    EXPECT_EQ(result.cols(), 3);

    // compute the true function values---should be exact (we are using 3rd order to estimate a degree 3 polynomial)
    const Eigen::Vector2d x_true = f(x);
    const Eigen::Vector2d y_true = f(y);
    const Eigen::Vector2d z_true = f(z);

    EXPECT_NEAR((x_true-result.col(0)).norm(), 0.0, 1.0e-10);
    EXPECT_NEAR((y_true-result.col(1)).norm(), 0.0, 1.0e-10);
    EXPECT_NEAR((z_true-result.col(2)).norm(), 0.0, 1.0e-10);
  }

  /// A matrix holding the input points.
  std::vector<Eigen::VectorXd> ins;

  /// A matrix holding the output points.
  std::vector<Eigen::VectorXd> outs;

  /// The order of the polynomial regression
  const unsigned int order = 3;

  /// The object that does the regression
  std::shared_ptr<Regression> reg;
};

TEST_F(RegressionTest, LegendreBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order);
}

TEST_F(RegressionTest, MonomialBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order, "Monomial");
}

TEST_F(RegressionTest, PhysicistHermiteBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order, "PhysicistHermite");
}

TEST_F(RegressionTest, ProbabilistHermiteBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order, "ProbabilistHermite");
}

TEST_F(RegressionTest, LaguerreBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order, "Laguerre");
}

TEST_F(RegressionTest, JacobiBasis) {
  // create the regression 
  reg = std::make_shared<Regression>(order, "Jacobi");
}

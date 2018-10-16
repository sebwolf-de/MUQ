#include <gtest/gtest.h>

#include "RosenbrockFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(CostFunctionTests, RosenbrockCost) {

  // the Rosenbrock cost function
  auto rosen = std::make_shared<RosenbrockFunction>();

  // choose a random point
  const Eigen::VectorXd x = Eigen::Vector2d::Random();
  const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 100.0);

  // the true value
  const double cst = (1.0-x(0))*(1.0-x(0))+100.0*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));

  // check the cost evaluations
  EXPECT_DOUBLE_EQ(cst, rosen->Evaluate(x, a).at(0) (0));
  EXPECT_DOUBLE_EQ(cst, rosen->Cost(ref_vector<Eigen::VectorXd>({std::cref(x), std::cref(a)})));
  EXPECT_DOUBLE_EQ(cst, rosen->Cost(x, a));

  // the true gradient
  const Eigen::Vector2d grad_true(-400.0*(x(1)-x(0)*x(0))*x(0)-2.0*(1.0-x(0)), 200.0*(x(1)-x(0)*x(0)));

  // compute the gradient
  const Eigen::VectorXd& grad_test0 = rosen->Gradient(0, x, a, (Eigen::VectorXd)Eigen::VectorXd::Ones(2));
  const Eigen::VectorXd& grad_test1 = rosen->Gradient(0, ref_vector<Eigen::VectorXd>({std::cref(x), std::cref(a)}), (Eigen::VectorXd)Eigen::VectorXd::Ones(2));

  EXPECT_DOUBLE_EQ((grad_true-grad_test0).norm(), 0.0);
  EXPECT_DOUBLE_EQ((grad_true-grad_test1).norm(), 0.0);


  // the true hessian
  Eigen::Matrix2d hess_temp;
  hess_temp << 1200.0*x(0)*x(0)-400.0*x(1)+2.0, -400.0*x(0),
               -400.0*x(0), 200.0;
  const Eigen::Matrix2d hess_true(hess_temp);

  // compute the Hessian
  std::vector<Eigen::VectorXd> input;
  input.push_back(x);
  input.push_back(a);

  const Eigen::MatrixXd& hess_test0 = 
    rosen->Hessian(0, input, (Eigen::VectorXd)Eigen::VectorXd::Ones(2));

  EXPECT_NEAR((hess_true-hess_test0).norm(), 0.0, 1.0e-5);

  // Test the Hessian action
  Eigen::Vector2d vec(-12.3, 34.6);


  const Eigen::VectorXd& hessAction_true = hess_true*vec;
  const Eigen::VectorXd& hessAction_test =
    rosen->ApplyHessian(0, input, (Eigen::VectorXd)Eigen::VectorXd::Ones(2), vec);
  
  EXPECT_NEAR((hessAction_true-hessAction_test).norm(), 0.0, 3.0e-4);
  
}


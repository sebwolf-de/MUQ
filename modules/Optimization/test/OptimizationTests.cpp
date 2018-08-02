#include <gtest/gtest.h>

#include "MUQ/Optimization/Optimization.h"

#include "RosenbrockFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(Optimization, BasicTest) {
  auto cost = std::make_shared<RosenbrockFunction>();

  // choose a random point
  const Eigen::VectorXd x = Eigen::Vector2d::Random();
  const boost::any xany = x;
  
  auto opt = std::make_shared<Optimization>(cost);

  std::vector<boost::any> soln0 = opt->Evaluate(x);
  std::pair<Eigen::VectorXd, double> soln1 = opt->Solve(x);
  std::pair<Eigen::VectorXd, double> soln2 = opt->Solve(ref_vector<boost::any>(1, std::cref(xany)));

  const Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd const&>(soln0[0]);
  const double minf = boost::any_cast<double const>(soln0[1]);

  EXPECT_EQ(xopt.size(), 2);
  EXPECT_NEAR(xopt(0), 1.0, 1.0e-5);
  EXPECT_NEAR(xopt(1), 1.0, 1.0e-5);
  EXPECT_NEAR(minf, 0.0, 1.0e-10);

  EXPECT_EQ(soln1.first.size(), 2);
  EXPECT_NEAR(soln1.first(0), 1.0, 1.0e-5);
  EXPECT_NEAR(soln1.first(1), 1.0, 1.0e-5);
  EXPECT_NEAR(soln1.second, 0.0, 1.0e-10);

  EXPECT_EQ(soln2.first.size(), 2);
  EXPECT_NEAR(soln2.first(0), 1.0, 1.0e-5);
  EXPECT_NEAR(soln2.first(1), 1.0, 1.0e-5);
  EXPECT_NEAR(soln2.second, 0.0, 1.0e-10);
}

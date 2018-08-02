#include <gtest/gtest.h>

#include "MUQ/Optimization/Optimization.h"

#include "RosenbrockFunction.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(Optimization, BasicTest) {
  pt::ptree pt;
  pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-14);
  pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-14);
  pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-14);
  pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-14);
  pt.put("Optimization.MaxEvaluations", 250); // max number of cost function evaluations
  
  auto cost = std::make_shared<RosenbrockFunction>();

  const Eigen::VectorXd x = Eigen::Vector2d(0.85, 1.2);
  const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);
  
  auto opt = std::make_shared<Optimization>(cost, pt.get_child("Optimization"));

  std::vector<boost::any> soln0 = opt->Evaluate(x, a);
  std::pair<Eigen::VectorXd, double> soln1 = opt->Solve(x, a);
  const boost::any xany = x;
  const boost::any aany = a;
  std::pair<Eigen::VectorXd, double> soln2 = opt->Solve(ref_vector<boost::any>({std::cref(xany), std::cref(aany)}));

  const Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd const&>(soln0[0]);
  const double minf = boost::any_cast<double const>(soln0[1]);

  EXPECT_EQ(xopt.size(), 2);
  EXPECT_NEAR(xopt(0), 1.0, 1.0e-4);
  EXPECT_NEAR(xopt(1), 1.0, 1.0e-4);
  EXPECT_NEAR(minf, 0.0, 1.0e-10);

  EXPECT_EQ(soln1.first.size(), 2);
  EXPECT_NEAR(soln1.first(0), 1.0, 1.0e-4);
  EXPECT_NEAR(soln1.first(1), 1.0, 1.0e-4);
  EXPECT_NEAR(soln1.second, 0.0, 1.0e-10);

  EXPECT_EQ(soln2.first.size(), 2);
  EXPECT_NEAR(soln2.first(0), 1.0, 1.0e-4);
  EXPECT_NEAR(soln2.first(1), 1.0, 1.0e-4);
  EXPECT_NEAR(soln2.second, 0.0, 1.0e-10);
}

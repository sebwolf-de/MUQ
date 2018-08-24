#include <gtest/gtest.h>

#include "MUQ/Optimization/Optimization.h"

#include "RosenbrockFunction.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;

class OptimizationTests : public::testing::Test {
public:

  inline OptimizationTests() {
    pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.MaxEvaluations", 1000); // max number of cost function evaluations
  }

  inline virtual ~OptimizationTests() {}
  
  inline void TearDown() {
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
  
  pt::ptree pt;
  
private:
};

TEST_F(OptimizationTests, COBYLA) {
  pt.put("Optimization.Algorithm", "COBYLA");
}

TEST_F(OptimizationTests, BOBYQA) {
  pt.put("Optimization.Algorithm", "BOBYQA");
}

TEST_F(OptimizationTests, NEWUOA) {
  pt.put("Optimization.Algorithm", "NEWUOA");
}

TEST_F(OptimizationTests, PRAXIS) {
  pt.put("Optimization.Algorithm", "PRAXIS");
}

TEST_F(OptimizationTests, NM) {
  pt.put("Optimization.Algorithm", "NM");
}

TEST_F(OptimizationTests, SBPLX) {
  pt.put("Optimization.Algorithm", "SBPLX");
}

TEST_F(OptimizationTests, MMA) {
  pt.put("Optimization.Algorithm", "MMA");
}

TEST_F(OptimizationTests, SLSQP) {
  pt.put("Optimization.Algorithm", "SLSQP");
}

TEST_F(OptimizationTests, LBFGS) {
  pt.put("Optimization.Algorithm", "LBFGS");
}

TEST_F(OptimizationTests, PreTN) {
  pt.put("Optimization.Algorithm", "PreTN");
}

TEST_F(OptimizationTests, LMVM) {
  pt.put("Optimization.Algorithm", "LMVM");
}

class InequalityConstraint : public CostFunction {
  public: 
  inline InequalityConstraint() : CostFunction(Eigen::Vector2i(2, 1)) {}

  virtual inline ~InequalityConstraint() {}

 private:

  inline virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {
    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd& b = input[1];

    return b(0)-xc(0);
  }

  inline virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override {
    assert(inputDimWrt==0);
    
    gradient = Eigen::Vector2d::Zero(2);
    gradient(0) = -1.0;

    gradient *= sensitivity(0);
  }
};

class EqualityConstraint : public CostFunction {
  public: 
  inline EqualityConstraint() : CostFunction(Eigen::VectorXi::Constant(1, 2)) {}

  virtual inline ~EqualityConstraint() {}

 private:

  inline virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {
    const Eigen::VectorXd& xc = input[0];

    return xc(1)-xc(0)*xc(0)-1.0;
  }

  inline virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override {
    assert(inputDimWrt==0);
    
    const Eigen::VectorXd& xc = input[0];

    gradient = Eigen::Vector2d::Zero(2);
    gradient(0) = -2.0*xc(0);
    gradient(1) = 1.0;

    gradient *= sensitivity(0);
  }
};

class ConstrainedOptimizationTests : public::testing::Test {
public:

  inline ConstrainedOptimizationTests() {
    pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.ConstraintTolerance", 1.0e-14);
    pt.put("Optimization.MaxEvaluations", 100000); // max number of cost function evaluations
  }

  inline virtual ~ConstrainedOptimizationTests() {}
  
  inline void TearDown() {
    auto cost = std::make_shared<RosenbrockFunction>();
    auto ineqconstraint = std::make_shared<InequalityConstraint>();
    auto eqconstraint = std::make_shared<EqualityConstraint>();
    
    const Eigen::VectorXd x = Eigen::Vector2d(2.0, 5.0);
    const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);
    const Eigen::VectorXd b = Eigen::VectorXd::Constant(1, 2.0);
	    
    auto opt = std::make_shared<Optimization>(cost, pt.get_child("Optimization"));
    opt->AddInequalityConstraint(ineqconstraint);
    opt->AddEqualityConstraint(eqconstraint);
    
    std::vector<boost::any> soln0 = opt->Evaluate(x, a, b);
    std::pair<Eigen::VectorXd, double> soln1 = opt->Solve(x, a, b);
    const boost::any xany = x;
    const boost::any aany = a;
    const boost::any bany = b;
    std::pair<Eigen::VectorXd, double> soln2 = opt->Solve(ref_vector<boost::any>({std::cref(xany), std::cref(aany), std::cref(bany)}));
    
    const Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd const&>(soln0[0]);
    const double minf = boost::any_cast<double const>(soln0[1]);

    EXPECT_EQ(xopt.size(), 2);
    EXPECT_NEAR(xopt(0), 2.0, 1.0e-4);
    EXPECT_NEAR(xopt(1), 5.0, 1.0e-4);
    EXPECT_NEAR(minf, 6.0, 1.0e-10);
    
    EXPECT_EQ(soln1.first.size(), 2);
    EXPECT_NEAR(soln1.first(0), 2.0, 1.0e-4);
    EXPECT_NEAR(soln1.first(1), 5.0, 1.0e-4);
    EXPECT_NEAR(soln1.second, 6.0, 1.0e-10);
    
    EXPECT_EQ(soln2.first.size(), 2);
    EXPECT_NEAR(soln2.first(0), 2.0, 1.0e-4);
    EXPECT_NEAR(soln2.first(1), 5.0, 1.0e-4);
    EXPECT_NEAR(soln2.second, 6.0, 1.0e-10);
  }
  
  pt::ptree pt;
  
private:
};

TEST_F(ConstrainedOptimizationTests, COBYLA) {
  pt.put("Optimization.Algorithm", "COBYLA");
}

TEST_F(ConstrainedOptimizationTests, SLSQP) {
  pt.put("Optimization.Algorithm", "SLSQP");
}

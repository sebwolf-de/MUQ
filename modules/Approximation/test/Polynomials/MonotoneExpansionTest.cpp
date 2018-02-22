#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "gtest/gtest.h"

using namespace muq::Approximation;
using namespace muq::Utilities;


class Approximation_MonotoneExpansion1d : public::testing::Test {
public:
  Approximation_MonotoneExpansion1d() {

    auto monomial = std::make_shared<Monomial>();
    auto bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(1, monomial);

    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1, 2);

    coeffs.resize(1,3);
    coeffs << -0.01, 0.1, -0.5; // intercept, slope

    monoPart = std::make_shared<BasisExpansion>(bases, multis, coeffs);
    expansion = std::make_shared<MonotoneExpansion>(monoPart);
  }

  virtual ~Approximation_MonotoneExpansion1d() {};

  std::shared_ptr<BasisExpansion> monoPart;
  std::shared_ptr<MonotoneExpansion> expansion;
  Eigen::MatrixXd coeffs;

};

TEST_F(Approximation_MonotoneExpansion1d, Evaluate){

  Eigen::VectorXd evalPt(1);

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt << lb;
  double oldOutput = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(evalPt).at(0))(0);

  for(int i=1; i<numSteps; ++i){
      evalPt(0) = lb + dx*double(i);
      double newOutput = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(evalPt).at(0))(0);
      EXPECT_GT(newOutput, oldOutput);

      double monoEval = boost::any_cast<Eigen::VectorXd>(monoPart->Evaluate(evalPt).at(0))(0);
      double trueMono = coeffs(0) + coeffs(1)*evalPt(0) + coeffs(2)*evalPt(0)*evalPt(0);
      EXPECT_DOUBLE_EQ(trueMono, monoEval);

      double truth = coeffs(0)*coeffs(0)*evalPt(0)
                     + coeffs(0)*coeffs(1)*std::pow(evalPt(0),2.0)
                     + (1.0/3.0)*(2.0*coeffs(0)*coeffs(2)+coeffs(1)*coeffs(1))*std::pow(evalPt(0),3.0)
                     + 0.5*coeffs(1)*coeffs(2)*std::pow(evalPt(0),4.0)
                     + 0.2*coeffs(2)*coeffs(2)*std::pow(evalPt(0),5.0);
      EXPECT_NEAR(truth, newOutput, 1e-2);

      oldOutput = newOutput;
  }
}

TEST_F(Approximation_MonotoneExpansion1d, EvaluateWithCoeffs){

  Eigen::VectorXd evalPt(1);

  Eigen::VectorXd newCoeffs(3);
  newCoeffs << 0.5, 2.0, -0.1;
  coeffs.row(0) = newCoeffs;

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt << lb;
  double oldOutput = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(evalPt, newCoeffs).at(0))(0);

  for(int i=1; i<numSteps; ++i){
      evalPt(0) = lb + dx*double(i);
      double newOutput = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(evalPt, newCoeffs).at(0))(0);
      EXPECT_GT(newOutput, oldOutput);

      double truth = coeffs(0)*coeffs(0)*evalPt(0)
                     + coeffs(0)*coeffs(1)*std::pow(evalPt(0),2.0)
                     + (1.0/3.0)*(2.0*coeffs(0)*coeffs(2)+coeffs(1)*coeffs(1))*std::pow(evalPt(0),3.0)
                     + 0.5*coeffs(1)*coeffs(2)*std::pow(evalPt(0),4.0)
                     + 0.2*coeffs(2)*coeffs(2)*std::pow(evalPt(0),5.0);
      EXPECT_NEAR(truth, newOutput, 5e-2);

      oldOutput = newOutput;
  }
}

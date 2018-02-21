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

    Eigen::MatrixXd coeffs(1,3);
    coeffs << 1.0, 2.0, 0.5; // intercept, slope

    expansion = std::make_shared<MontoneExpansion>(std::make_shared<BasisExpansion>(bases, multis, coeffs));
  }

  virtual ~Approximation_MonotoneExpansion1d() {};

  std::shared_ptr<MonotoneExpansion> expansion;

};

TEST_F(Approximation_MonotoneExpansion1d, Evaluate){

  Eigen::VectorXd evalPt(1);
  evalPt << 0.5;

  Eigen::VectorXd output = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(evalPt)[0]);

  std::cout << "Output = " << output(0) << std::endl;

}

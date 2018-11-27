#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/OptimalExperimentalDesign/Evidence.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

class Likelihood : public Distribution {
public:

  inline Likelihood() : Distribution(1, Eigen::VectorXi::Ones(2)) {
    gauss = std::make_shared<Gaussian>(1);
  }

  virtual ~Likelihood() = default;
private:

  virtual inline double LogDensityImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    return gauss->LogDensity(inputs);
  }

  virtual inline Eigen::VectorXd SampleImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    return gauss->Sample(inputs);
  }

  std::shared_ptr<Gaussian> gauss;
};

TEST(EvidenceTest, Basic) {
  auto prior = std::make_shared<Gaussian>(1);
  auto like = std::make_shared<Likelihood>();

  pt::ptree pt;
  pt.put("NumImportanceSamples", 200);
  auto evidence = std::make_shared<Evidence>(prior, like, prior, pt);

  const Eigen::VectorXd y = Eigen::VectorXd::Random(1); // data
  const Eigen::VectorXd d = Eigen::VectorXd::Random(1); // design

  const double piy = evidence->LogDensity(y, d);

  auto expectedEvidence = std::make_shared<Gaussian>(1);

  const double piyExpected = expectedEvidence->LogDensity(y);

  EXPECT_NEAR(piy, piyExpected, 1.0e-10);
}

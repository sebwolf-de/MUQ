#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/OptimalExperimentalDesign/Evidence.h"
#include "MUQ/OptimalExperimentalDesign/OEDResidual.h"

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

  virtual inline double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    return gauss->LogDensity(inputs);
  }

  virtual inline Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    return gauss->Sample(inputs);
  }

  std::shared_ptr<Gaussian> gauss;
};

TEST(OEDResidualTest, Basic) {
  auto prior = std::make_shared<Gaussian>(1);
  auto like = std::make_shared<Likelihood>();

  pt::ptree pt;
  pt.put("Evidence.NumImportanceSamples", 200);
  auto evidence = std::make_shared<Evidence>(prior, like, prior, pt.get_child("Evidence"));

  pt.put("OEDResidual.NumImportanceSamples", 200);
  auto resid = std::make_shared<OEDResidual>(like, evidence, like, pt.get_child("OEDResidual"));

  const Eigen::VectorXd x = Eigen::VectorXd::Random(1); // parameter
  const Eigen::VectorXd d = Eigen::VectorXd::Random(1); // design

  const Eigen::VectorXd result = resid->Evaluate(x, d) [0];
  EXPECT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0), 0.0, 1.0e-10);
}

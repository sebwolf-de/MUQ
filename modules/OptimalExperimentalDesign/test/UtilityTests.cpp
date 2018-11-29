#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/OptimalExperimentalDesign/Evidence.h"
#include "MUQ/OptimalExperimentalDesign/OEDResidual.h"
#include "MUQ/OptimalExperimentalDesign/Utility.h"

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

class UtilityBiasing : public Distribution {
public:
  UtilityBiasing(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood) : Distribution(2, Eigen::VectorXi::Ones(1)), prior(prior), likelihood(likelihood) {}

  virtual ~UtilityBiasing() = default;
private:
  virtual inline Eigen::VectorXd SampleImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    const Eigen::VectorXd& d = inputs[0];
    const Eigen::VectorXd& x = prior->Sample();
    const Eigen::VectorXd& y = likelihood->Sample(x, d);

    return Eigen::Vector2d(x(0), y(0));
  }


  virtual inline double LogDensityImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    const Eigen::VectorXd& x = inputs[0].get().head(1);
    const Eigen::VectorXd& y = inputs[0].get().tail(1);
    const Eigen::VectorXd& d = inputs[1];

    return prior->LogDensity(x)+likelihood->LogDensity(y, x, d);
  }

  std::shared_ptr<muq::Modeling::Distribution> prior;

  std::shared_ptr<muq::Modeling::Distribution> likelihood;
};

TEST(UtilityTest, Basic) {
  auto prior = std::make_shared<Gaussian>(1);
  auto like = std::make_shared<Likelihood>();

  pt::ptree pt;
  pt.put("Evidence.NumImportanceSamples", 200);
  auto evidence = std::make_shared<Evidence>(prior, like, prior, pt.get_child("Evidence"));

  auto biasing = std::make_shared<UtilityBiasing>(prior, like);

  pt.put("Utility.NumImportanceSamples", 200);
  auto utility = std::make_shared<Utility>(prior, like, evidence, biasing, pt.get_child("Utility"));

  const Eigen::VectorXd d = Eigen::VectorXd::Random(1); // design

  const Eigen::VectorXd result = utility->Evaluate(d) [0];
  EXPECT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0), 0.0, 1.0e-10);
}

TEST(UtilityTest, Surrogate) {
  auto prior = std::make_shared<Gaussian>(1);
  auto like = std::make_shared<Likelihood>();

  pt::ptree pt;
  pt.put("Evidence.NumImportanceSamples", 25);
  auto evidence = std::make_shared<Evidence>(prior, like, prior, pt.get_child("Evidence"));

  pt.put("OEDResidual.NumImportanceSamples", 25);
  auto resid = std::make_shared<OEDResidual>(like, evidence, like, pt.get_child("OEDResidual"));

  pt.put("Utility.NumImportanceSamples", 25);
  pt.put("Utility.LocalRegression", "MyLocalReg");
  pt.put("Utility.MyLocalReg.NumNeighbors", 11);
  pt.put("Utility.MyLocalReg.Order", 2);
  auto utility = std::make_shared<Utility>(prior, resid, prior, pt.get_child("Utility"));

  const Eigen::VectorXd d = Eigen::VectorXd::Random(1); // design

  const Eigen::VectorXd result = utility->Evaluate(d) [0];
  EXPECT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0), 0.0, 1.0e-10);
}

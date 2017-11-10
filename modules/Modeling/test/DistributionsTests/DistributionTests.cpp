#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Distribution.h"

using namespace muq::Modeling;

class ExampleDensity : public Distribution {
public:

  inline ExampleDensity() : Distribution() {}

  inline virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) const override {
    // get the point where we are evaluating the log density
    const double x = boost::any_cast<double const>(inputs[0]);

    return -x*x-std::sin(2.0*M_PI*x);
  }
  
private:
};

TEST(Distribution, EvaluateDensity) {
  // create a distribution that only has a (log) density
  auto dens = std::make_shared<ExampleDensity>();

  // evaluate the density at a point
  const double x = 2.5;

  // evaluate the log density
  double logdens = dens->LogDensity(x);

  // make sure we get the density we implemented
  EXPECT_DOUBLE_EQ(logdens, -x*x-std::sin(2.0*M_PI*x));

  // evaluate the log density using evaluate
  const auto& result = dens->Evaluate(Distribution::Mode::EvaluateLogDensity, x);
  logdens = boost::any_cast<double const>(result[0]);

  // make sure we get the density we implemented
  EXPECT_DOUBLE_EQ(logdens, -x*x-std::sin(2.0*M_PI*x));
}

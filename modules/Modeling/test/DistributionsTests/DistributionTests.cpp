#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

class ExampleDensity : public Distribution {
public:

  inline ExampleDensity() : Distribution() {}

  inline virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) override {
    // get the point where we are evaluating the log density
    const double x = boost::any_cast<double const>(inputs[0]);

    return -x*x-std::sin(2.0*M_PI*x);
  }

private:
};

class ExampleRV : public Distribution {
public:
    inline ExampleRV() : Distribution() {}

    inline virtual boost::any SampleImpl(ref_vector<boost::any> const& inputs) override {
        double outputVal = 0.1;
        return boost::any(outputVal);
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
  std::vector<boost::any> const& result = dens->Evaluate(Distribution::Mode::EvaluateLogDensity, x);
  logdens = boost::any_cast<double const>(result[0]);
  // make sure we get the density we implemented
  EXPECT_DOUBLE_EQ(logdens, -x*x-std::sin(2.0*M_PI*x));

  std::shared_ptr<Density> densPiece = dens->AsDensity();
  ASSERT_TRUE(densPiece);

  std::vector<boost::any> const& result2 = densPiece->Evaluate(x);
  double logDens2 = boost::any_cast<double const>(result.at(0));
  EXPECT_DOUBLE_EQ(logDens2, -x*x-std::sin(2.0*M_PI*x));

  // Make sure we can't sample
  EXPECT_THROW(dens->Sample(x), muq::NotImplementedError);
  EXPECT_THROW(dens->Evaluate(Distribution::Mode::SampleDistribution, x), muq::NotImplementedError);
}

TEST(Distribution, EvaluateSample) {

    auto rv = std::make_shared<ExampleRV>();
    const double x  = 2.5;

    // draw a sample
    boost::any samp = rv->Sample(x);
    EXPECT_DOUBLE_EQ(0.1, boost::any_cast<double>(samp));

    samp = rv->Evaluate(Distribution::Mode::SampleDistribution, x)[0];
    EXPECT_DOUBLE_EQ(0.1, boost::any_cast<double>(samp));

    // Make sure we can't evaluate the density
    EXPECT_THROW(rv->LogDensity(x), muq::NotImplementedError);
    EXPECT_THROW(rv->Evaluate(Distribution::Mode::EvaluateLogDensity, x), muq::NotImplementedError);

}

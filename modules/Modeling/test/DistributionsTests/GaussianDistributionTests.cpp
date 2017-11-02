#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

using namespace muq::Modeling;

class GaussianDistributionTests : public::testing::Test {
public:

  inline GaussianDistributionTests() {}

  inline virtual ~GaussianDistributionTests() {}
  
private:
};

TEST_F(GaussianDistributionTests, EvaluateLogDensity) {
  EXPECT_TRUE(false);
}

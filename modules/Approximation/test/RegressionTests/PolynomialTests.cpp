#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/Monomial.h"
#include "MUQ/Approximation/Regression/Hermite.h"
#include "MUQ/Approximation/Regression/Legendre.h"

using namespace muq::Approximation;

TEST(Polynomial, Monomial) {
  // create a monomial object
  auto mono = std::make_shared<Monomial>();
  EXPECT_EQ(mono->numInputs, 2);
  EXPECT_EQ(mono->numOutputs, 1);

  // a point to evaluate the monomial
  const double x = 2.0;

  // the order of the polynomial
  const unsigned int max_order = 200;

  for( unsigned int p=0; p<max_order; ++p ) {
    // evaluate the monomial
    const std::vector<boost::any>& result = mono->Evaluate(p, x);
    EXPECT_EQ(result.size(), 1);
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(result[0]), std::pow(x, p));
  }
}

TEST(Polynomial, Hermite) {
  // create a Hermite object
  auto hermite = std::make_shared<Hermite>();
  
  EXPECT_DOUBLE_EQ(1.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)0, 0.4) [0]));
  EXPECT_DOUBLE_EQ(0.6, boost::any_cast<double const>(hermite->Evaluate((unsigned int)1, 0.3) [0]));
  EXPECT_DOUBLE_EQ(4.0*std::pow(0.6, 2.0)-2.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)2, 0.6) [0]));
  EXPECT_DOUBLE_EQ(33.6235290625, boost::any_cast<double const>(hermite->Evaluate((unsigned int)5, 0.325) [0]));
  EXPECT_NEAR(6219.5581337600015, boost::any_cast<double const>(hermite->Evaluate((unsigned int)8, 1.6) [0]), 2.0e-12);
  EXPECT_NEAR(6.075804453410837e11, boost::any_cast<double const>(hermite->Evaluate((unsigned int)20, -0.845) [0]), 3.0e-4);
}

TEST(Polynomial, Legendre) {
  // create a Legendre object
  auto legendre = std::make_shared<Legendre>();

  // A sample of points tested against mathematica
  EXPECT_DOUBLE_EQ(1.0, boost::any_cast<double const>(legendre->Evaluate((unsigned int)0, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.3, boost::any_cast<double const>(legendre->Evaluate((unsigned int)1, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.5*(3.0*std::pow(0.3, 2.0)-1.0), boost::any_cast<double const>(legendre->Evaluate((unsigned int)2, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.3375579333496094, boost::any_cast<double const>(legendre->Evaluate((unsigned int)5, 0.325) [0]));
  EXPECT_NEAR(-0.05346106275520913, boost::any_cast<double const>(legendre->Evaluate((unsigned int)20, -0.845) [0]), 1.0e-14);
  EXPECT_NEAR(-0.1119514835092105, boost::any_cast<double const>(legendre->Evaluate((unsigned int)50, 0.1264) [0]), 1e-14);
  EXPECT_NEAR(-0.001892916076323403, boost::any_cast<double const>(legendre->Evaluate((unsigned int)200, -0.3598) [0]), 1e-14);
  EXPECT_NEAR(0.01954143166718206, boost::any_cast<double const>(legendre->Evaluate((unsigned int)1000, 0.4587) [0]), 1e-14);
}

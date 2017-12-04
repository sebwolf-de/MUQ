#include <gtest/gtest.h>

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

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

    // Get the first derivative
    double deriv = mono->DerivativeEvaluate(p,1,x);
    EXPECT_DOUBLE_EQ(double(p)*std::pow(x,p-1.0), deriv);

    // Get the second derivative
    deriv = mono->DerivativeEvaluate(p,2,x);
    EXPECT_DOUBLE_EQ(double(p)*double(p-1)*std::pow(x,p-2.0), deriv);

    // Third derivative
    deriv = mono->DerivativeEvaluate(p,3,x);
    EXPECT_DOUBLE_EQ(double(p)*double(p-1)*double(p-2)*std::pow(x,p-3.0), deriv);
    
  }
}

TEST(Polynomial, PhysicistHermite) {
  // create a Hermite object
  auto hermite = std::make_shared<PhysicistHermite>();

  // Evaluations
  EXPECT_DOUBLE_EQ(1.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)0, 0.4) [0]));
  EXPECT_DOUBLE_EQ(0.6, boost::any_cast<double const>(hermite->Evaluate((unsigned int)1, 0.3) [0]));
  EXPECT_DOUBLE_EQ(4.0*std::pow(0.6, 2.0)-2.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)2, 0.6) [0]));
  EXPECT_DOUBLE_EQ(33.6235290625, boost::any_cast<double const>(hermite->Evaluate((unsigned int)5, 0.325) [0]));
  EXPECT_NEAR(6219.5581337600015, boost::any_cast<double const>(hermite->Evaluate((unsigned int)8, 1.6) [0]), 2.0e-12);
  EXPECT_NEAR(6.075804453410837e11, boost::any_cast<double const>(hermite->Evaluate((unsigned int)20, -0.845) [0]), 3.0e-4);

  // First derivatives
  const double x = 0.23;
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(2.0, hermite->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(8.0*x, hermite->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(24.0*x*x - 12.0, hermite->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(64.0*std::pow(x,3.0) -96.0*x, hermite->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(8.0, hermite->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(48.0*x, hermite->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(192.0*std::pow(x,2.0) - 96.0, hermite->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(48.0, hermite->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(384.0*x, hermite->DerivativeEvaluate(4,3,x));  
  
}

TEST(Polynomial, ProbabilistHermite) {
  // create a Hermite object
  auto hermite = std::make_shared<ProbabilistHermite>();

  const double x = 1.32;
  EXPECT_DOUBLE_EQ(1.0, hermite->PolynomialEvaluate(0, x));
  EXPECT_DOUBLE_EQ(x, hermite->PolynomialEvaluate(1, x));
  EXPECT_DOUBLE_EQ(x*x-1.0, hermite->PolynomialEvaluate(2, x));
  EXPECT_DOUBLE_EQ(x*x*x - 3.0*x, hermite->PolynomialEvaluate(3, x));
  EXPECT_DOUBLE_EQ(x*x*x*x - 6.0*x*x + 3.0, hermite->PolynomialEvaluate(4, x));
  EXPECT_NEAR(std::pow(x,5) - 10*std::pow(x,3) + 15.0*x, hermite->PolynomialEvaluate(5, x), 1e-10);

  // First derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(1.0, hermite->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(2.0*x, hermite->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(3.0*x*x - 3.0, hermite->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(4.0*std::pow(x,3.0) - 12.0*x, hermite->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(2.0, hermite->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(6.0*x, hermite->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(12.0*std::pow(x,2.0) - 12.0, hermite->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(6.0, hermite->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(24.0*x, hermite->DerivativeEvaluate(4,3,x)); 

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

  // First derivatives
  const double x = 0.23;
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(1.0, legendre->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(3.0*x, legendre->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(7.5*x*x - 1.5, legendre->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(17.5*std::pow(x,3.0) - 7.5*x, legendre->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(3.0, legendre->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(15.0*x, legendre->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(52.5*std::pow(x,2.0) - 7.5, legendre->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(15.0, legendre->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(105.0*x, legendre->DerivativeEvaluate(4,3,x));
}

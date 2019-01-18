#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

using namespace muq::Approximation;

///Regression test for a 6-dim mixed var type and mixed order full tensor quadrature grid
TEST(UtilitiesFullTensorQuadrature, regressionTest)
{
  auto hermitePoly = std::make_shared<PhysicistHermite>();
  auto gaussHermite = std::make_shared<GaussQuadrature>(hermitePoly,5);

  auto legendrePoly = std::make_shared<PhysicistHermite>();
  auto gaussLegendre = std::make_shared<GaussQuadrature>(legendrePoly,5);

  FullTensorQuadrature quad({gaussHermite, gaussLegendre});
  EXPECT_EQ(gaussLegendre->Points().cols()*gaussHermite->Points().cols(), quad.Points().cols());

  EXPECT_NEAR(3.14159, quad.Weights().sum(),2e-5);

}
